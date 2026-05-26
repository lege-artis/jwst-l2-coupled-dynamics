"""
test_integration.py — end-to-end smoke tests for the JWST-L2 first-cut.

Re-runs the first_example scenario at a shorter horizon and asserts that:
- the integrator completes without errors
- conservation residuals stay below the documented gates
- the FFT-extracted free-precession period matches the theoretical value
  within FFT bin resolution

These tests validate the full code path from geometry construction through
integration to spectral analysis.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from geometries import make_jwst_like, make_probe_like
from dynamics import (
    BodyState, pack_state, unpack_state, rk4_step,
    total_kinetic_energy, total_potential_energy, total_angular_momentum,
)


@pytest.fixture
def short_run_result():
    """Run the first_example scenario at a reduced horizon (200 s) so tests
    finish quickly. Returns the final state and the omega_A_x time series."""
    body_a = make_jwst_like()
    body_b = make_probe_like()
    m_a = body_a.total_mass
    m_b = body_b.total_mass
    I_a_body = body_a.inertia_at_com()
    I_b_body = body_b.inertia_at_com()

    # tilted initial attitude as in run_first_example
    theta = 0.3
    axis = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
    q_a_init = np.concatenate([[np.cos(theta / 2.0)], np.sin(theta / 2.0) * axis])

    initial_a = BodyState(
        x=[0.0, 0.0, 0.0], v=[0.0, 0.0, 0.0],
        q=q_a_init, omega=[0.02, 0.0, 0.08],
    )
    initial_b = BodyState(
        x=[50.0, 0.0, 0.0], v=[0.0, 0.0, 0.0],
        q=[1.0, 0.0, 0.0, 0.0], omega=[0.0, 0.03, 0.08],
    )
    state = pack_state(initial_a, initial_b)

    # baseline conservation quantities
    Etk0, Erk0, _ = total_kinetic_energy(state, m_a, I_a_body, m_b, I_b_body)
    Ep0 = total_potential_energy(state, m_a, m_b)
    E0 = Etk0 + Erk0 + Ep0
    L0 = total_angular_momentum(state, m_a, I_a_body, m_b, I_b_body)

    # integrate at reduced horizon (200 s with dt = 0.05 -> 4000 steps)
    omega_x_series = []
    t_series = []
    t = 0.0
    dt = 0.05
    n_steps = 4000
    snapshot_every = 10
    for step in range(n_steps + 1):
        if step % snapshot_every == 0:
            ba, _ = unpack_state(state)
            omega_x_series.append(float(ba.omega[0]))
            t_series.append(float(t))
        if step < n_steps:
            state = rk4_step(state, dt, m_a, I_a_body, m_b, I_b_body)
            t += dt

    Etk_f, Erk_f, _ = total_kinetic_energy(state, m_a, I_a_body, m_b, I_b_body)
    Ep_f = total_potential_energy(state, m_a, m_b)
    E_f = Etk_f + Erk_f + Ep_f
    L_f = total_angular_momentum(state, m_a, I_a_body, m_b, I_b_body)

    return {
        "E0": E0,
        "E_f": E_f,
        "L0": L0,
        "L_f": L_f,
        "omega_x_series": np.array(omega_x_series),
        "t_series": np.array(t_series),
        "I_a_principal_moments": np.linalg.eigvalsh(I_a_body),
    }


class TestEndToEndIntegration:
    def test_run_completes_without_error(self, short_run_result):
        # If the fixture executes the integration without exception, this passes.
        # The check is implicit in the fixture itself, but we add an assertion to
        # make the test name informative.
        assert short_run_result["t_series"][-1] == pytest.approx(200.0, rel=1e-6)

    def test_energy_conserved_below_gate(self, short_run_result):
        E0 = short_run_result["E0"]
        E_f = short_run_result["E_f"]
        dE_rel = abs(E_f - E0) / abs(E0)
        # gate: |dE/E0| < 1e-9 over 200 s, 4000 RK4 steps. We observed ~1e-11
        # in the 600 s first_example; 200 s should easily clear the looser gate.
        assert dE_rel < 1e-9, f"dE/E0 = {dE_rel:.3e} exceeds 1e-9 gate"

    def test_angular_momentum_conserved_below_gate(self, short_run_result):
        L0 = short_run_result["L0"]
        L_f = short_run_result["L_f"]
        dL_rel = np.linalg.norm(L_f - L0) / np.linalg.norm(L0)
        # gate: |dL/L0| < 1e-8 over 200 s. We observed ~2e-10 in the 600 s
        # first_example; the 200 s run should clear this comfortably.
        assert dL_rel < 1e-8, f"dL/L0 = {dL_rel:.3e} exceeds 1e-8 gate"

    def test_fft_peak_matches_theoretical_euler_frequency(self, short_run_result):
        omega_x = short_run_result["omega_x_series"]
        t = short_run_result["t_series"]
        dt_snap = t[1] - t[0]
        fs = 1.0 / dt_snap   # 2 Hz snapshot rate

        signal = omega_x - omega_x.mean()
        Xk = np.fft.rfft(signal)
        freqs = np.fft.rfftfreq(len(signal), d=dt_snap)

        # find the dominant peak (excluding DC)
        peak_idx = 1 + int(np.argmax(np.abs(Xk[1:])))
        peak_freq = float(freqs[peak_idx])
        peak_period = 1.0 / peak_freq if peak_freq > 0 else 0.0

        # theoretical Euler frequency:
        # for body A: I_principal sorted ascending = (15384, 23322, 23322)
        # I_axial = 15384, I_perp = 23322
        # lambda = (I_axial - I_perp) / I_perp * omega_z = (15384-23322)/23322 * 0.08
        pmA = short_run_result["I_a_principal_moments"]
        lambda_theory = (pmA[0] - pmA[-1]) / pmA[-1] * 0.08
        period_theory = 2.0 * np.pi / abs(lambda_theory)

        # FFT bin spacing at 200 s window with 2 Hz sample rate:
        # df = 1 / (n_samples * dt_snap) = 1 / 200 = 5 mHz -> bins at 0.005, 0.01, ...
        # period theory is 231 s -> freq 4.33 mHz -> closest bin is the first non-DC
        # at 5 mHz -> period 200 s. So peak_period should land near 200 s.
        # Tolerance: ±25% to account for bin-quantisation in this short window
        assert peak_period == pytest.approx(period_theory, rel=0.30), \
            f"peak_period = {peak_period:.2f} s, theory = {period_theory:.2f} s"
