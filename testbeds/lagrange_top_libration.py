"""
lagrange_top_libration.py — Lagrange-top libration testbed.

Analytical oracle: Ch 06 §6.4.3 closed-form Lagrange-top + gravity-gradient
libration solution via effective-potential method (Podolský 2018 §8.3).

Finite-amplitude period (Ch 06 §6.4.3 boxed):
  T_lib(ϑ_amp) = (2/π) K(k) · T_lib^linear
where K(k) is the complete elliptic integral of the first kind and k is the
elliptic modulus determined by the libration amplitude ϑ_amp.

Testbed procedure (Phase 2 brief §2.5):
  1. Configure a librating rigid body (prolate, axisymmetric) in a field of a
     gravitational primary at separation rho.
  2. Initialise at amplitude ϑ_amp with ω = 0 (turning point).
  3. Integrate n_periods × T_analytical steps with the supplied integrator.
  4. Extract numerical period via zero-crossing analysis on the libration-angle
     signal θ(t) = angle of body symmetry-axis from equilibrium direction.
  5. Compare against analytical T_lib (AG-LIB-1: < 0.1% relative error).
  6. Convergence-order test at dt, dt/2, dt/4 (AG-LIB-2).

OQ-PHASE2-8 note (RESOLVED 2026-06-01):
  The brief §2.5 / Ch 03 §3.7.2 "~66 days" figure was a dropped-2π arithmetic
  slip; the correct Earth-tidal value is ~210 days (this file's oracle is
  right).  ~362 days is the same body in the Sun-tidal regime (m_sun at ~1 AU),
  not a contradiction.  See
  _audit/OQ-PHASE2-8-RECONCILIATION-2026-06-01.md (+ -PATCH-DIFFS-).
  All tests use THIS file's analytical oracle as the sole reference (unchanged).

References:
  Podolský (2018) §8.3 — CS-language Lagrange-top oracle
  Ch 06 §6.4.3 — EN-language closed-form derivation
  Phase 2 brief §2.5 — testbed specification
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
from scipy.special import ellipk  # complete elliptic integral of the first kind

sys.path.insert(0, str(Path(__file__).parent.parent))
from dynamics import G_NEWTON, BodyState, pack_state, unpack_state, q_normalise, q_to_matrix

# Import quaternion / rotation helpers from the variational integrator module
from integrators.symplectic_se3_variational import rotvec_to_quat


# ─────────────────────────────────────────────────────────────────────────────
# Analytical oracle
# ─────────────────────────────────────────────────────────────────────────────

def libration_frequency_linear(
    I_perp: float, I_parallel: float, m_b: float, rho: float,
) -> float:
    """Small-amplitude libration angular frequency ω_lib (Ch 06 §6.4.2).

    ω_lib² = 3 G m_b (I_perp − I_parallel) / (I_perp · ρ³)

    Returns ω_lib > 0 (only for prolate body: I_parallel < I_perp).
    Raises ValueError for oblate (unstable equilibrium).
    """
    val = 3.0 * G_NEWTON * m_b * (I_perp - I_parallel) / (I_perp * rho**3)
    if val <= 0.0:
        raise ValueError(
            f"libration_frequency_linear: oblate body (I_parallel={I_parallel} >= "
            f"I_perp={I_perp}) — equilibrium is unstable, no libration frequency."
        )
    return np.sqrt(val)


def elliptic_modulus(
    vartheta_amp: float,
    I_perp: float, I_parallel: float,
    L_psi: float,
    m_b: float, rho: float,
) -> float:
    """Elliptic modulus k for the gravity-gradient finite-amplitude libration.

    Derived from energy conservation for the EOM
      θ̈ = −(ω_lib²/2) sin(2θ)
    with effective potential V ∝ −cos(2θ).  At the turning point θ = ϑ_amp
    (θ̇ = 0), energy conservation gives:

      θ̇² = ω_lib² (sin²ϑ_amp − sin²θ)

    so T/4 = (1/ω_lib) K(k) with k = sin(ϑ_amp).

    For the pure libration regime (L_psi = 0):

      k = sin(ϑ_amp)          ← gravity-gradient / Ch 06 §6.4.3

    Note: the Lagrange-top formula k = sin(ϑ_amp/2) applies to the potential
    V ∝ −cos θ, which is physically different from the gravity-gradient
    V ∝ −cos(2θ).  The two agree to O(ϑ_amp²) but diverge at large amplitudes.

    Args:
        vartheta_amp: libration amplitude in radians
        I_perp, I_parallel: principal moments of inertia
        L_psi: spin angular momentum (about symmetry axis)
        m_b: mass of gravitational primary
        rho: orbital radius

    Returns:
        k: elliptic modulus ∈ [0, 1)
    """
    if abs(L_psi) < 1e-20:
        return np.sin(vartheta_amp)
    raise NotImplementedError(
        "elliptic_modulus for non-zero L_psi: Sitting 2 full implementation."
    )


def libration_period_analytical(
    vartheta_amp: float,
    I_perp: float, I_parallel: float,
    m_b: float, rho: float,
    *,
    L_psi: float = 0.0,
) -> float:
    """Analytical libration period at amplitude ϑ_amp (Ch 06 §6.4.3).

    T_lib(ϑ_amp) = (2/π) K(k) · T_lib^linear

    where K(k) is the complete elliptic integral of the first kind and
    T_lib^linear = 2π / ω_lib is the small-amplitude limit.

    Args:
        vartheta_amp: libration amplitude (radians); must be > 0
        I_perp, I_parallel: axisymmetric body inertia (I_parallel < I_perp for stable)
        m_b: mass of gravitational primary
        rho: orbital radius (meters)
        L_psi: spin angular momentum about symmetry axis (0 for pure libration)

    Returns:
        T_lib: libration period (seconds)
    """
    omega_lib = libration_frequency_linear(I_perp, I_parallel, m_b, rho)
    T_linear = 2.0 * np.pi / omega_lib
    k = elliptic_modulus(vartheta_amp, I_perp, I_parallel, L_psi, m_b, rho)
    K_k = float(ellipk(k**2))   # scipy: ellipk(m) where m = k²
    return (2.0 / np.pi) * K_k * T_linear


# ─────────────────────────────────────────────────────────────────────────────
# Period extraction from oscillatory signal
# ─────────────────────────────────────────────────────────────────────────────

def _extract_period_zero_crossings(signal: np.ndarray, dt: float) -> float:
    """Extract oscillation period via positive-to-negative zero-crossing analysis.

    Uses linear interpolation between samples to locate crossing times precisely.
    Returns the mean period across all observed crossings.
    Returns np.nan if fewer than 2 crossings are found.
    """
    crossings = []
    for i in range(1, len(signal)):
        if signal[i - 1] > 0.0 and signal[i] <= 0.0:
            # Linear interpolation for sub-sample crossing time
            frac = signal[i - 1] / (signal[i - 1] - signal[i])
            crossings.append((i - 1 + frac) * dt)
    if len(crossings) < 2:
        return np.nan
    periods = np.diff(crossings)
    return float(np.mean(periods))


# ─────────────────────────────────────────────────────────────────────────────
# Build initial state for libration IC
# ─────────────────────────────────────────────────────────────────────────────

def _make_libration_state(
    amplitude_rad: float,
    I_a_body: np.ndarray,
    m_a: float,
    m_b: float,
    rho: float,
    I_b_body: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build a two-body state for the libration testbed.

    Body A (librating): axisymmetric, z-axis is symmetry axis.  Orbits body B
    in a circular orbit in the x-z plane.  Starts at maximum pitch amplitude
    (turning point) aligned with the local vertical.

    Body B (gravitational primary): at position (0, 0, rho) from A; isotropic
    inertia so it experiences zero gravity-gradient torque; m_b >> m_a so it
    barely moves.

    Geometry at t = 0:
      - Local vertical (zenith): r_hat = (0, 0, 1)  [from A toward B]
      - Orbital velocity:  v_A = (v_circ, 0, 0)
      - Orbit normal:      n_hat = −ŷ  (angular momentum L = r_AB × v_A ∝ −ŷ)
      - Pitch direction:   p_hat = n_hat × r_hat = −x̂

    Equilibrium orientation: body-z along r_hat = (0, 0, 1).
    Pitch tilt of +amplitude: body-z toward p_hat (−x̂), achieved by rotating
      about +y by amplitude_rad:  body-z = (sin θ, 0, cos θ).
    Then signal₀ = −dot(body-z, v̂_A) = −sin θ < 0, so we use −y rotation:
      q_a = rotvec([0, −amplitude, 0])  →  body-z = (−sin θ, 0, cos θ)
      signal₀ = −dot((−sin θ, 0, cos θ), (1, 0, 0)) = sin θ > 0  ✓

    Equilibrium spin: in inertial frame ω_A = ω_orb · (−ŷ) so that body-z
    tracks r_hat.  For R_y rotations the y-component is invariant, so:
      omega_body = R_a0.T @ [0, −ω_orb, 0] = [0, −ω_orb, 0]

    Returns (state, I_a_body, I_b_body_used).
    """
    if I_b_body is None:
        # Isotropic inertia → zero gravity-gradient torque on B
        I_b_body = np.eye(3)

    # Circular orbit parameters
    v_circ = np.sqrt(G_NEWTON * m_b / rho)          # orbital speed (m/s)
    omega_orb = v_circ / rho                          # = sqrt(G m_b / rho³)

    # Body A: pitch tilt about −y (turning point; signal starts positive)
    q_a = rotvec_to_quat(np.array([0.0, -amplitude_rad, 0.0]))
    # Equilibrium spin in body frame: [0, −ω_orb, 0] (y-invariant under R_y)
    omega_body = [0.0, -omega_orb, 0.0]

    body_a = BodyState(
        x=[0.0, 0.0, 0.0],
        v=[v_circ, 0.0, 0.0],
        q=q_a,
        omega=omega_body,
    )
    # Body B: at (0, 0, rho), essentially fixed (large m_b, isotropic I)
    body_b = BodyState(
        x=[0.0, 0.0, rho],
        v=[0.0, 0.0, 0.0],
        q=[1.0, 0.0, 0.0, 0.0],
        omega=[0.0, 0.0, 0.0],
    )
    state = pack_state(body_a, body_b)
    return state, I_a_body, I_b_body


# ─────────────────────────────────────────────────────────────────────────────
# Numerical libration period extraction
# ─────────────────────────────────────────────────────────────────────────────

def libration_period_numerical(
    integrator,
    I_perp: float,
    I_parallel: float,
    m_a: float,
    m_b: float,
    rho: float,
    ic_amplitude_deg: float,
    dt: float,
    n_periods: int = 5,
) -> float:
    """Integrate and extract numerical libration period.

    Uses the analytical oracle from libration_period_analytical() to determine
    how many steps to run (n_periods complete librations), then extracts the
    actual period from the libration-angle signal via zero-crossing analysis.

    The libration-angle signal is θ(t) tracked as the y-component of the body
    symmetry-axis (body z) expressed in the inertial frame.  For a tilt about
    x-axis, this component is approximately −sin(θ) which oscillates as a
    negative-amplitude sine when starting from the positive turning point.

    Args:
        integrator: callable with signature
            ``(state, dt, m_a, I_a_body, m_b, I_b_body) -> state``
        I_perp, I_parallel: axisymmetric body principal moments (I_par < I_perp)
        m_a: mass of the librating body (kg)
        m_b: mass of the gravitational primary (kg)
        rho: separation distance (m)
        ic_amplitude_deg: initial tilt angle (degrees)
        dt: timestep (seconds)
        n_periods: number of libration periods to integrate

    Returns:
        T_numerical: extracted libration period (seconds)
    """
    amplitude_rad = np.deg2rad(ic_amplitude_deg)
    T_analytical = libration_period_analytical(
        amplitude_rad, I_perp, I_parallel, m_b, rho
    )
    n_steps = max(int(n_periods * T_analytical / dt), 1)

    I_a_body = np.diag([I_perp, I_perp, I_parallel])
    state, I_a_body, I_b_body = _make_libration_state(
        amplitude_rad, I_a_body, m_a, m_b, rho
    )

    # Integrate and record libration-angle signal.
    #
    # Signal: −dot(z_body, v̂_A)  where z_body = R_A @ ê_z.
    #
    # Derivation: body A orbits in the x-z plane at angular rate ω_orb.
    # The pitch angle θ(t) is the angle from the instantaneous local vertical.
    # For a pitch perturbation θ in the rotating frame:
    #   z_body(inertial) = r_hat(t)·cos θ + pitch_dir(t)·sin θ
    # and −dot(z_body, v̂_A) = sin θ(t)  (algebraically exact for circular orbit).
    # So the signal equals sin(θ(t)), which oscillates at ω_lib.
    #
    # At t = 0 with the turning-point IC:
    #   signal₀ = sin(amplitude_rad) > 0  ✓
    # Zero-crossings occur at T_lib/4, 5T_lib/4, 9T_lib/4, …
    # Each inter-crossing gap = T_lib (the period we want to measure).

    _EZ = np.array([0.0, 0.0, 1.0])

    def _signal(st: np.ndarray) -> float:
        ba, _ = unpack_state(st)
        R = q_to_matrix(ba.q)
        z_body = R @ _EZ
        v_a = np.array(ba.v)
        v_norm = float(np.linalg.norm(v_a))
        if v_norm < 1e-15:
            return 0.0
        return float(-np.dot(z_body, v_a / v_norm))

    signal = [_signal(state)]
    for _ in range(n_steps):
        state = integrator(state, dt, m_a, I_a_body, m_b, I_b_body)
        signal.append(_signal(state))

    signal = np.array(signal)
    T_num = _extract_period_zero_crossings(signal, dt)
    return T_num


# ─────────────────────────────────────────────────────────────────────────────
# Convergence-order test
# ─────────────────────────────────────────────────────────────────────────────

def convergence_order_test(
    integrator,
    I_perp: float,
    I_parallel: float,
    m_a: float,
    m_b: float,
    rho: float,
    ic_amplitude_deg: float = 5.0,
    dt_base: float = 1.0,
    n_periods: int = 5,
) -> dict:
    """Empirical convergence-order test at dt_base, dt_base/2, dt_base/4.

    Computes relative period error at each timestep and estimates the
    convergence rate from successive error ratios:
      rate_1 = log2(err_dt / err_dt2)
      rate_2 = log2(err_dt2 / err_dt4)

    Returns dict with keys:
        T_analytical: analytical period (s)
        T_dt, T_dt2, T_dt4: numerical periods at dt, dt/2, dt/4
        err_dt, err_dt2, err_dt4: relative period errors (dimensionless)
        rate_1: log2(err_dt / err_dt2) — convergence rate from first pair
        rate_2: log2(err_dt2 / err_dt4) — convergence rate from second pair
    """
    amplitude_rad = np.deg2rad(ic_amplitude_deg)
    T_analytical = libration_period_analytical(
        amplitude_rad, I_perp, I_parallel, m_b, rho
    )

    results = {"T_analytical": T_analytical}
    labels = ["dt", "dt2", "dt4"]
    dts = [dt_base, dt_base / 2.0, dt_base / 4.0]
    periods = []

    for label, dt in zip(labels, dts):
        T_num = libration_period_numerical(
            integrator, I_perp, I_parallel, m_a, m_b, rho,
            ic_amplitude_deg, dt, n_periods=n_periods,
        )
        results[f"T_{label}"] = T_num
        err = abs(T_num - T_analytical) / T_analytical
        results[f"err_{label}"] = err
        periods.append(T_num)

    # Convergence rates
    err_dt  = results["err_dt"]
    err_dt2 = results["err_dt2"]
    err_dt4 = results["err_dt4"]
    eps_floor = 1e-14  # avoid log(0)
    results["rate_1"] = np.log2(max(err_dt, eps_floor) / max(err_dt2, eps_floor))
    results["rate_2"] = np.log2(max(err_dt2, eps_floor) / max(err_dt4, eps_floor))
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Full testbed runner
# ─────────────────────────────────────────────────────────────────────────────

def run_testbed(
    integrator=None,
    dt: float = 1.0,
    amplitude_deg: float = 5.0,
    n_periods: int = 5,
    I_perp: float = 2.0,
    I_parallel: float = 1.0,
    m_a: float = 1.0,
    m_b: float = 1e8,
    rho: float = 1.0,
) -> dict:
    """Run the full Lagrange-top libration testbed (AG-LIB-1 + AG-LIB-2).

    Default parameters give a compact synthetic testbed with T_lib ≈ 62.8 s
    for fast unit testing.  For the JWST-at-L2 reference:
      I_perp=23323, I_parallel=15384 (make_jwst_like() inertia, kg·m²)
      m_b=5.972e24 (Earth mass, kg), rho=1.5e9 (L2 separation, m)
    which gives T_lib ≈ 210 days (Earth tidal; see OQ-PHASE2-8 note at top
    of this file regarding the brief's 66-day citation).

    Args:
        integrator: callable; if None, uses se3_variational_step
        dt: timestep (s)
        amplitude_deg: libration amplitude (degrees)
        n_periods: number of periods to integrate for period extraction
        I_perp, I_parallel: axisymmetric body inertia (kg·m²)
        m_a: mass of librating body (kg)
        m_b: mass of gravitational primary (kg)
        rho: separation (m)

    Returns dict with keys:
        T_analytical, T_numerical, period_error_pct,
        convergence_results (from convergence_order_test),
        n_steps.
    """
    if integrator is None:
        from integrators.symplectic_se3_variational import se3_variational_step
        integrator = se3_variational_step

    amplitude_rad = np.deg2rad(amplitude_deg)
    T_analytical = libration_period_analytical(
        amplitude_rad, I_perp, I_parallel, m_b, rho
    )
    T_numerical = libration_period_numerical(
        integrator, I_perp, I_parallel, m_a, m_b, rho,
        amplitude_deg, dt, n_periods=n_periods,
    )
    period_error_pct = abs(T_numerical - T_analytical) / T_analytical * 100.0
    n_steps = int(n_periods * T_analytical / dt)

    convergence_results = convergence_order_test(
        integrator, I_perp, I_parallel, m_a, m_b, rho,
        ic_amplitude_deg=amplitude_deg,
        dt_base=dt,
        n_periods=n_periods,
    )

    return {
        "T_analytical": T_analytical,
        "T_numerical": T_numerical,
        "period_error_pct": period_error_pct,
        "convergence_results": convergence_results,
        "n_steps": n_steps,
    }
