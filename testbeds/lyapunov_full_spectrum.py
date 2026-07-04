"""
lyapunov_full_spectrum.py — 12-vector full-spectrum testbed (Phase 3).

Runs the full 12-exponent Lyapunov spectrum on the three calibration scenarios
using both:
  Tier A — quaternion-deflated R^14 cocycle (lyapunov.cocycle_deflated)
  Tier B — Lie-algebraic so(3)×R^3 cocycle  (lyapunov.cocycle_liealg)

Scenario initial conditions are identical to testbeds/lyapunov_three_scenario.py
(DO NOT EDIT that file); this testbed is additive.

Production run parameters for a definitive Phase 3 measurement:
  Tier A scenario (c):  dt=50 s, n_steps=4000  → T=200 000 s  (≈28 flip periods)
  Tier B scenario (a):  dt=1 s,  n_steps=20000 → T=20 000 s   (≈220 flip periods)
  Tier B scenario (c):  dt=50 s, n_steps=2000  → T=100 000 s  (≈14 flip periods)

Fast-run defaults (for CI / structure checks) use much smaller n_steps; the
production_run=True flag switches to the larger parameters.

PHASE 3 DELIVERABLE.

References:
  testbeds/lyapunov_three_scenario.py   (scenario IC — DO NOT EDIT)
  lyapunov/cocycle_deflated.py          (Tier A)
  lyapunov/cocycle_liealg.py            (Tier B)
  lyapunov/pairing.py                   (gate helpers)
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

# ── Shared constants (mirror lyapunov_three_scenario.py) ─────────────────────
_M_EARTH = 5.97e24
_R_EARTH = 6.371e6
_G = 6.67430e-11


# ─────────────────────────────────────────────────────────────────────────────
# Scenario-setup helpers (reuse IC logic from lyapunov_three_scenario.py)
# ─────────────────────────────────────────────────────────────────────────────

def _scenario_a_state_and_params():
    """IC + params for scenario (a): torque-free asymmetric top."""
    from geometries import make_probe_like
    from dynamics import BodyState, pack_state

    I_a = np.diag([100.0, 200.0, 400.0])
    m_a = 500.0
    omega0 = np.array([0.1, 0.0, 0.0])  # stable axis

    probe = make_probe_like()
    m_b = probe.total_mass
    I_b = probe.inertia_at_com()

    body_a = BodyState(np.zeros(3), np.zeros(3),
                       np.array([1.0, 0.0, 0.0, 0.0]), omega0)
    body_b = BodyState(np.array([1e12, 0.0, 0.0]), np.zeros(3),
                       np.array([1.0, 0.0, 0.0, 0.0]), np.zeros(3))

    return pack_state(body_a, body_b), m_a, I_a, m_b, I_b


def _scenario_c_state_and_params():
    """IC + params for scenario (c): GG-Dzhanibekov near-saddle."""
    from dynamics import BodyState, pack_state

    I1, I2, I3 = 100.0, 200.0, 400.0
    I_a = np.diag([I1, I2, I3])
    m_a = 500.0
    omega0 = np.array([1.0e-3, 1.0e-2, 1.0e-3])

    m_b = _M_EARTH
    I_b = np.eye(3) * (2.0 / 5.0 * m_b * _R_EARTH**2)
    rho = 3.0e6

    body_a = BodyState(np.zeros(3), np.zeros(3),
                       np.array([1.0, 0.0, 0.0, 0.0]), omega0)
    body_b = BodyState(np.array([rho, 0.0, 0.0]), np.zeros(3),
                       np.array([1.0, 0.0, 0.0, 0.0]), np.zeros(3))

    return pack_state(body_a, body_b), m_a, I_a, m_b, I_b


# ─────────────────────────────────────────────────────────────────────────────
# Tier A — deflated R^14 spectrum
# ─────────────────────────────────────────────────────────────────────────────

def run_tier_a_scenario_c(
    dt: float = 50.0,
    n_steps: int = 4000,
    n_vectors: int = 12,
    *,
    production_run: bool = False,
) -> dict:
    """Deflated 12-vector spectrum, scenario (c).  AG-P3-1 + AG-P3-2 target.

    Production run (production_run=True): dt=50 s, n_steps=4000 → T=200 000 s.
    Fast run default: use caller-supplied dt / n_steps.
    """
    from integrators.symplectic_se3_variational import se3_variational_step
    from lyapunov.cocycle_deflated import lyapunov_spectrum_deflated
    from lyapunov.pairing import full_spectrum_gates
    from lyapunov.reference import dzhanibekov_characteristic_frequency

    if production_run:
        dt, n_steps = 50.0, 4000

    state0, m_a, I_a, m_b, I_b = _scenario_c_state_and_params()

    result = lyapunov_spectrum_deflated(
        state0, se3_variational_step, dt, n_steps, n_vectors,
        m_a, I_a, m_b, I_b,
    )

    gates = full_spectrum_gates(result['exponents'])
    sigma = dzhanibekov_characteristic_frequency(100.0, 200.0, 400.0, 1e-2)

    return {
        'tier': 'A',
        'scenario': 'c',
        'exponents': result['exponents'],
        'sum_exponents': result['sum_exponents'],
        'sigma': sigma,
        'gates': gates,
        'T_total': result['T_total'],
        'n_qr': result['n_qr'],
        'tau_qr_used': result['tau_qr_used'],
        'running_exponents': result.get('running_exponents', []),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Tier B — Lie-algebraic spectrum
# ─────────────────────────────────────────────────────────────────────────────

def run_tier_b_scenario_a(
    dt: float = 1.0,
    n_steps: int = 20_000,
    n_vectors: int = 1,
    *,
    production_run: bool = False,
) -> dict:
    """Lie-alg λ_max, scenario (a): AG-P3-3 target.

    Production run: dt=1 s, n_steps=20 000 → T=20 000 s.
    """
    from integrators.symplectic_se3_variational import se3_variational_step
    from lyapunov.cocycle_liealg import lyapunov_spectrum_liealg
    from lyapunov.reference import lambda_max_gate_scenario_a

    if production_run:
        dt, n_steps = 1.0, 20_000

    state0, m_a, I_a, m_b, I_b = _scenario_a_state_and_params()

    result = lyapunov_spectrum_liealg(
        state0, se3_variational_step, dt, n_steps, n_vectors,
        m_a, I_a, m_b, I_b,
    )

    lambda_max = float(result['exponents'][0])
    # AG-P3-3 aspirational gate: < 1e-5 (vs 2.3e-4 in R^14).
    # Document the empirical value; gate will be confirmed post-measurement.
    ag_p3_3_aspirational = lambda_max < 1e-5
    ag_p3_3_phase2_parity = lambda_max_gate_scenario_a(lambda_max)  # < 1e-3

    return {
        'tier': 'B',
        'scenario': 'a',
        'exponents': result['exponents'],
        'lambda_max': lambda_max,
        'ag_p3_3_aspirational_pass': ag_p3_3_aspirational,
        'ag_p3_3_phase2_parity_pass': ag_p3_3_phase2_parity,
        'T_total': result['T_total'],
        'n_qr': result['n_qr'],
        'tau_qr_used': result['tau_qr_used'],
        'running_exponents': result.get('running_exponents', []),
    }


def run_tier_b_scenario_c(
    dt: float = 50.0,
    n_steps: int = 2000,
    n_vectors: int = 12,
    *,
    production_run: bool = False,
) -> dict:
    """Lie-alg full spectrum, scenario (c): AG-P3-2 + AG-P3-4 target.

    Production run: dt=50 s, n_steps=2 000 → T=100 000 s.
    """
    from integrators.symplectic_se3_variational import se3_variational_step
    from lyapunov.cocycle_liealg import lyapunov_spectrum_liealg
    from lyapunov.pairing import full_spectrum_gates
    from lyapunov.reference import dzhanibekov_characteristic_frequency

    if production_run:
        dt, n_steps = 50.0, 2000

    state0, m_a, I_a, m_b, I_b = _scenario_c_state_and_params()

    result = lyapunov_spectrum_liealg(
        state0, se3_variational_step, dt, n_steps, n_vectors,
        m_a, I_a, m_b, I_b,
    )

    gates = full_spectrum_gates(result['exponents'])
    sigma = dzhanibekov_characteristic_frequency(100.0, 200.0, 400.0, 1e-2)

    return {
        'tier': 'B',
        'scenario': 'c',
        'exponents': result['exponents'],
        'sum_exponents': result['sum_exponents'],
        'sigma': sigma,
        'gates': gates,
        'T_total': result['T_total'],
        'n_qr': result['n_qr'],
        'tau_qr_used': result['tau_qr_used'],
        'running_exponents': result.get('running_exponents', []),
    }
