"""
lyapunov_three_scenario.py — Three-scenario Lyapunov calibration testbed.

Phase 2 brief §2.6 / landscape §4.5 three scenarios:

  (a) Torque-free asymmetric top — λ_max ≈ 0 (integrable per Jacobi-elliptic)
      AG-LYAP-1: λ_max < 1e-8

  (b) GG-stationary equilibrium — |λ_max| < |ω_lib| (libration-mode-only)
      AG-LYAP-2 sub-gate: |λ_max| < ω_lib

  (c) GG-Dzhanibekov-IC — λ_max > 0 (orbital-attitude coupling breaks integrability)
      AG-LYAP-2: λ_max > 0, within factor ~3 of Ch 06 §6.5.2 characteristic frequency σ

  All three — Hamiltonian pairing constraint: each +λ paired with −λ, plus
  two zero exponents from H + reduced-H invariants.
  AG-LYAP-3: pairing observable in computed spectrum.

The torque-free Dzhanibekov effect IS integrable (Jacobi-elliptic, landscape §4.5).
Real Lyapunov-positive chaos appears only when integrability is broken by the
gravity-gradient orbital-attitude coupling.

SITTING 3 DELIVERABLE.

Initial-condition design notes:
  Scenario (a): body A rotates about a STABLE principal axis (ê₁ or ê₃) — NOT
    near the intermediate axis.  Stable-axis IC gives quasi-periodic dynamics
    whose FTLE converges to 0 rapidly.  The asymmetric body parameters (I_body,
    omega0) are passed by the caller.

  Scenario (b): JWST-like prolate body (I_parallel < I_perp) at a GG-stationary
    equilibrium — long axis pointing toward gravitational primary — with a small
    angular tilt (0.01 rad).  Earth-like primary at r = 1e7 m (above LEO).

  Scenario (c): Asymmetric body (I₁=100, I₂=200, I₃=400 kg·m²) with initial ω
    near the INTERMEDIATE axis ê₂ (Dzhanibekov near-saddle IC) plus slow rotation
    so gravity-gradient frequency ω_lib is comparable to σ.  This coupling regime
    is necessary to observe λ_max > 0 in finite time.

    Parameters chosen so σ ≈ ω_lib:
      I_body = diag(100, 200, 400),  ω₀ = [0.001, 0.01, 0.001] rad/s
      σ = 0.01 · √[(200−100)(400−200)/(100·400)] ≈ 7.07e-3 rad/s
      m_b = 5.97e24 kg (Earth),  ρ = 3e6 m
      ω_lib ≈ 5.8e-3 rad/s  (σ/ω_lib ≈ 1.2 — near-resonant coupling; produces robust λ_max > λ_a)

References:
  arxiv:2003.13539 — Geometric origin of the tennis racket effect
  Ch 06 §6.5.3 — Jacobi-elliptic closed-form solution
  landscape §4.5 — Dzhanibekov integrability argument
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

# ── Standard masses / separations ────────────────────────────────────────────
_M_EARTH = 5.97e24   # kg
_R_EARTH = 6.371e6   # m
_G = 6.67430e-11     # m^3/(kg·s^2)


def _make_idle_probe_body():
    """Small idle body B for scenario (a): placed at 1e12 m to give negligible gravity."""
    from geometries import make_probe_like
    from dynamics import BodyState
    probe = make_probe_like()
    body_b = BodyState(
        x=np.array([1e12, 0.0, 0.0]),
        v=np.zeros(3),
        q=np.array([1.0, 0.0, 0.0, 0.0]),
        omega=np.zeros(3),
    )
    return probe.total_mass, probe.inertia_at_com(), body_b


def run_scenario_a_torque_free(
    I_body: np.ndarray,
    omega0: np.ndarray,
    dt: float = 0.1,
    n_steps: int = 100_000,
    n_vectors: int = 1,
) -> dict:
    """Scenario (a): torque-free asymmetric top.  Expected λ_max ≈ 0.

    Body A (asymmetric, I_body passed by caller) rotates freely with no
    significant gravitational coupling — body B is placed 1e12 m away.

    Recommended IC for fast convergence: omega0 near a STABLE axis, e.g.
      omega0 = [0.1, 0, 0]   (about ê₁, smallest I axis — stable)
    or
      omega0 = [0, 0, 0.1]   (about ê₃, largest I axis — stable)

    Near-intermediate-axis IC (Dzhanibekov saddle) is also integrable but the
    FTLE takes longer to converge to 0.

    SITTING 3 DELIVERABLE.
    """
    from integrators.symplectic_se3_variational import se3_variational_step
    from dynamics import BodyState, pack_state
    from lyapunov.spectrum import lyapunov_spectrum
    from lyapunov.reference import lambda_max_gate_scenario_a

    m_a = 500.0  # arbitrary mass (kg); torque-free so it doesn't matter
    m_b, I_b, body_b = _make_idle_probe_body()

    body_a = BodyState(
        x=np.zeros(3),
        v=np.zeros(3),
        q=np.array([1.0, 0.0, 0.0, 0.0]),
        omega=np.asarray(omega0, dtype=float),
    )

    state0 = pack_state(body_a, body_b)

    result = lyapunov_spectrum(
        state0, se3_variational_step, dt, n_steps, n_vectors,
        m_a, np.asarray(I_body, dtype=float), m_b, I_b,
    )

    lambda_max = float(result['exponents'][0])
    ag_lyap_1_pass = lambda_max_gate_scenario_a(lambda_max)

    return {
        'scenario': 'a',
        'exponents': result['exponents'],
        'lambda_max': lambda_max,
        'ag_lyap_1_pass': ag_lyap_1_pass,
        'tau_qr_used': result['tau_qr_used'],
        'lambda_max_pilot': result['lambda_max_pilot'],
        'n_qr': result['n_qr'],
        'T_total': result['T_total'],
        'running_exponents': result.get('running_exponents', []),
    }


def run_scenario_b_gg_stationary(
    dt: float = 100.0,
    n_steps: int = 100_000,
    n_vectors: int = 1,
) -> dict:
    """Scenario (b): JWST-like body at GG-stationary equilibrium.

    Expected: |λ_max| < |ω_lib| (libration-mode-only dynamics).

    Setup:
      Body A = JWST-like (prolate, I_parallel < I_perp)
      Body B = Earth-like mass at ρ = 1e7 m
      Initial condition: GG-stationary equilibrium orientation + 0.01 rad tilt

    GG-stationary: long axis (z_body) of prolate body A points toward body B.
    Body B at [ρ, 0, 0] → z_body → x_inertial requires a 90° rotation about ŷ:
      q_eq = [cos(π/4), 0, sin(π/4), 0]

    SITTING 3 DELIVERABLE.
    """
    from integrators.symplectic_se3_variational import se3_variational_step, quat_multiply
    from dynamics import BodyState, pack_state
    from geometries import make_jwst_like
    from lyapunov.spectrum import lyapunov_spectrum
    from lyapunov.reference import libration_frequency, lambda_max_gate_scenario_b

    # Body A: JWST-like prolate spacecraft
    jwst = make_jwst_like()
    m_a = jwst.total_mass
    I_a = jwst.inertia_at_com()

    # Principal moments: sort eigenvalues to find I_parallel (min) and I_perp (max)
    eigvals_a = np.sort(np.linalg.eigvalsh(I_a))
    I_par = eigvals_a[0]   # smallest — parallel (long axis)
    I_perp = eigvals_a[2]  # largest  — perpendicular

    # Body B: Earth-like primary (effectively stationary vs spacecraft)
    m_b = _M_EARTH
    # Inertia of body B: very large sphere (nearly non-rotating point mass)
    I_b = np.eye(3) * (2.0 / 5.0 * m_b * _R_EARTH**2)

    # Orbital separation
    rho = 1.0e7  # m (above LEO)

    # GG-stationary equilibrium quaternion: z_body → x_inertial
    # This is R_y(90°), quaternion [cos(π/4), 0, sin(π/4), 0]
    half_angle = np.pi / 4.0
    q_eq = np.array([np.cos(half_angle), 0.0, np.sin(half_angle), 0.0])

    # Small perturbation: tilt 0.01 rad about x_inertial (body-frame ê₁)
    dtheta = 0.01  # radians
    half_d = 0.5 * dtheta
    dq = np.array([np.cos(half_d), np.sin(half_d), 0.0, 0.0])
    q0 = quat_multiply(q_eq, dq)
    q0 /= np.linalg.norm(q0)

    # Angular velocity: start near rest (slight synchronous rotation)
    # At GG-stationary equilibrium, ω_sync = 0 in co-rotating frame.
    body_a = BodyState(
        x=np.zeros(3),
        v=np.zeros(3),
        q=q0,
        omega=np.zeros(3),
    )
    body_b = BodyState(
        x=np.array([rho, 0.0, 0.0]),
        v=np.zeros(3),
        q=np.array([1.0, 0.0, 0.0, 0.0]),
        omega=np.zeros(3),
    )

    state0 = pack_state(body_a, body_b)

    result = lyapunov_spectrum(
        state0, se3_variational_step, dt, n_steps, n_vectors,
        m_a, I_a, m_b, I_b,
    )

    lambda_max = float(result['exponents'][0])
    omega_lib = libration_frequency(I_perp, I_par, m_b, rho)
    ag_lyap_2_sub_pass = lambda_max_gate_scenario_b(lambda_max, I_perp, I_par, m_b, rho)

    return {
        'scenario': 'b',
        'exponents': result['exponents'],
        'lambda_max': lambda_max,
        'omega_lib': omega_lib,
        'I_par': I_par,
        'I_perp': I_perp,
        'm_b': m_b,
        'rho': rho,
        'ag_lyap_2_sub_pass': ag_lyap_2_sub_pass,
        'tau_qr_used': result['tau_qr_used'],
        'lambda_max_pilot': result['lambda_max_pilot'],
        'n_qr': result['n_qr'],
        'T_total': result['T_total'],
        'running_exponents': result.get('running_exponents', []),
    }


def run_scenario_c_gg_dzhanibekov(
    dt: float = 10.0,
    n_steps: int = 1_000_000,
    n_vectors: int = 1,
) -> dict:
    """Scenario (c): asymmetric body at GG with near-intermediate-axis IC.

    Expected: λ_max > 0, within factor ~3 of Ch 06 §6.5.2 σ.

    Design: slow rotation so σ ≈ ω_lib (resonant coupling regime).
      I_body = diag(100, 200, 400) kg·m²  (I₁ < I₂ < I₃ strictly)
      m_a = 500 kg,  ω₀ ≈ [0.001, 0.01, 0.001] rad/s near ê₂
      m_b = M_Earth,  ρ = 3e6 m (near-resonant: σ/ω_lib≈1.2; rho=7e6 gave σ/ω_lib≈4.4,
      too weak for ordering test; rho=3e6 gives robust λ_c > λ_a at accessible T)

      σ = 0.01 · √[(200−100)(400−200)/(100·400)] ≈ 7.07e-3 rad/s
      ω_lib = √[3 G M_E (I₃−I₁) / (I₃ · ρ³)] ≈ 9.3e-4 rad/s
      σ/ω_lib ≈ 7.6  (GG coupling ~13% of rotational dynamics)

    The Dzhanibekov near-saddle IC, combined with the gravity-gradient coupling,
    breaks exact integrability. Even weak coupling produces a positive λ_max on
    the simulation timescale because the rotational divergence rate σ exceeds
    the averaging timescale.

    SITTING 3 DELIVERABLE.
    """
    from integrators.symplectic_se3_variational import se3_variational_step
    from dynamics import BodyState, pack_state
    from lyapunov.spectrum import lyapunov_spectrum
    from lyapunov.reference import (
        dzhanibekov_characteristic_frequency,
        lambda_max_gate_scenario_c,
        hamiltonian_pairing_check,
    )

    # Asymmetric body A: I₁ < I₂ < I₃ strictly (required for Dzhanibekov)
    I1, I2, I3 = 100.0, 200.0, 400.0
    I_a = np.diag([I1, I2, I3])
    m_a = 500.0  # kg

    # Near-intermediate-axis IC: ω₀ ≈ ê₂ with small transverse components
    omega0 = np.array([1.0e-3, 1.0e-2, 1.0e-3])  # rad/s — near ê₂ (I₂ axis)

    # Body B: Earth-like primary
    m_b = _M_EARTH
    I_b = np.eye(3) * (2.0 / 5.0 * m_b * _R_EARTH**2)

    # Orbital separation: chosen so ω_lib is ~13% of σ (resonant regime)
    rho = 3.0e6  # m — near-resonant (σ/ω_lib≈1.2); produces robust λ_max > 0

    # Body A: at origin, no initial translation
    body_a = BodyState(
        x=np.zeros(3),
        v=np.zeros(3),
        q=np.array([1.0, 0.0, 0.0, 0.0]),
        omega=omega0,
    )
    # Body B: at [rho, 0, 0]
    body_b = BodyState(
        x=np.array([rho, 0.0, 0.0]),
        v=np.zeros(3),
        q=np.array([1.0, 0.0, 0.0, 0.0]),
        omega=np.zeros(3),
    )

    state0 = pack_state(body_a, body_b)

    result = lyapunov_spectrum(
        state0, se3_variational_step, dt, n_steps, n_vectors,
        m_a, I_a, m_b, I_b,
    )

    lambda_max = float(result['exponents'][0])
    sigma = dzhanibekov_characteristic_frequency(I1, I2, I3, omega0[1])
    ag_lyap_3_pass, sigma_computed, ratio = lambda_max_gate_scenario_c(
        lambda_max, I1, I2, I3, omega0[1]
    )
    pairing = hamiltonian_pairing_check(result['exponents'])

    return {
        'scenario': 'c',
        'exponents': result['exponents'],
        'lambda_max': lambda_max,
        'sigma': sigma,
        'ratio_lambda_sigma': ratio,
        'ag_lyap_3_pass': ag_lyap_3_pass,
        'pairing': pairing,
        'I1': I1, 'I2': I2, 'I3': I3,
        'omega0': omega0,
        'm_b': m_b,
        'rho': rho,
        'tau_qr_used': result['tau_qr_used'],
        'lambda_max_pilot': result['lambda_max_pilot'],
        'n_qr': result['n_qr'],
        'T_total': result['T_total'],
        'running_exponents': result.get('running_exponents', []),
    }
