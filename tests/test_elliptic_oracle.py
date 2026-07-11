"""
test_elliptic_oracle.py — regression guard for Ch 06 §6.5.3 (corrected 2026-07-11).

Validates the DOCUMENTED boxed Jacobi-elliptic parametrisation of the torque-free
asymmetric top against a tight-tolerance solve_ivp integration of the exact Euler
equations. This is the cross-check that TC-DZH-02's fixture metadata claimed but
never implemented — its absence let erroneous boxed formulas reach the public
mirror (v0.1.0..v0.3.0). See _audit/forensic-2026-07-11/ and the §6.5.3
correction note.

Guard contract: the formulas below are transcribed from the CORRECTED Ch 06
§6.5.3 boxes. If anyone edits the doc boxes, they must edit these in lockstep —
a mismatch against numerics fails here first (KB-043 discipline: tests are
physical-truth-anchored, not documentation-anchored).

Branch covered: |L|^2 > 2*T*I2 (polhode encircles e3 — the §6.5.5 Dzhanibekov
test configurations). Runtime: ~10 s.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.special import ellipk, ellipj, ellipkinc

I1, I2, I3 = 100.0, 200.0, 400.0


def euler_rhs(t, w):
    wx, wy, wz = w
    return [(I2 - I3) * wy * wz / I1,
            (I3 - I1) * wz * wx / I2,
            (I1 - I2) * wx * wy / I3]


def invariants(w0):
    L = np.array([I1, I2, I3]) * np.asarray(w0)
    return float(L @ L), float(I1 * w0[0] ** 2 + I2 * w0[1] ** 2 + I3 * w0[2] ** 2)


def ch06_closed_form(M2, twoE):
    """Corrected Ch 06 §6.5.3 boxed forms (branch |L|^2 > 2T·I2)."""
    assert M2 > twoE * I2, "config must lie on the near-e3 branch"
    A = np.sqrt(I1 * (twoE * I3 - M2) / (I3 - I1))
    B = np.sqrt(I2 * (twoE * I3 - M2) / (I3 - I2))
    C = np.sqrt(I3 * (M2 - twoE * I1) / (I3 - I1))
    Omega_star = np.sqrt((I3 - I2) * (M2 - twoE * I1) / (I1 * I2 * I3))
    m = (I2 - I1) * (twoE * I3 - M2) / ((I3 - I2) * (M2 - twoE * I1))  # m = k^2
    T_LL = 4.0 * ellipk(m) / Omega_star
    return A, B, C, Omega_star, m, T_LL


def golden(w0, t_max, n_samp=200_000):
    sol = solve_ivp(euler_rhs, (0.0, t_max), w0, rtol=1e-12, atol=1e-14,
                    dense_output=True, max_step=t_max / 1e5)
    ts = np.linspace(0.0, t_max, n_samp)
    return ts, sol.sol(ts)


def flip_interval(ts, wy):
    idx = np.where(np.diff(np.sign(wy)) != 0)[0]
    cross = [ts[i] - wy[i] * (ts[i + 1] - ts[i]) / (wy[i + 1] - wy[i]) for i in idx]
    assert len(cross) >= 3, "need >=3 zero crossings in window"
    return float(np.mean(np.diff(cross)))


W0_MAIN = np.array([1e-3, 1e-2, 1e-3])          # TC-DZH-02 configuration


class TestEllipticOracle:
    def test_flip_period_closed_form_vs_golden(self):
        """T_LL/2 (documented closed form) == mean omega_2 zero-crossing interval."""
        M2, twoE = invariants(W0_MAIN)
        *_, T_LL = ch06_closed_form(M2, twoE)
        ts, w = golden(W0_MAIN, t_max=8000.0)
        meas = flip_interval(ts, w[1])
        assert abs(T_LL / 2.0 - meas) / meas < 1e-9

    def test_amplitude_is_turning_point_B(self):
        """max |L_y| equals the documented amplitude B (the old boxed alpha was ~16x off)."""
        M2, twoE = invariants(W0_MAIN)
        _, B, *_ = ch06_closed_form(M2, twoE)
        ts, w = golden(W0_MAIN, t_max=4000.0)
        Ly_max = float(np.max(np.abs(I2 * w[1])))
        assert abs(Ly_max - B) / B < 1e-6

    def test_parametrisation_elementwise(self):
        """(L_x, L_y, L_z)(t) from the documented boxes matches solve_ivp element-wise."""
        M2, twoE = invariants(W0_MAIN)
        A, B, C, Om, m, T_LL = ch06_closed_form(M2, twoE)
        t_max = 2.0 * T_LL
        ts, w = golden(W0_MAIN, t_max=t_max, n_samp=20_000)
        Lnum = np.array([I1, I2, I3])[:, None] * w
        # phase from initial L_y via F(arcsin(sn), m); resolve branch by best fit
        s0 = np.clip(I2 * W0_MAIN[1] / B, -1.0, 1.0)
        u0 = ellipkinc(np.arcsin(s0), m)
        K = ellipk(m)
        best = np.inf
        for uc in (u0, 2 * K - u0, -u0, -(2 * K - u0)):
            for sgn in (+1.0, -1.0):
                sn, cn, dn, _ = ellipj(Om * ts + uc, m)
                Lpar = np.vstack([sgn * A * cn, sgn * B * sn, C * dn])
                best = min(best, np.max(np.abs(Lpar - Lnum)) / np.max(np.abs(Lnum)))
        assert best < 1e-6

    def test_conservation_identities_exact(self):
        """A^2+C^2 relations: parametrisation satisfies both conservation laws identically."""
        M2, twoE = invariants(W0_MAIN)
        A, B, C, Om, m, _ = ch06_closed_form(M2, twoE)
        u = np.linspace(0.0, 4 * ellipk(m), 1000)
        sn, cn, dn, _ = ellipj(u, m)
        L = np.vstack([A * cn, B * sn, C * dn])
        M2_t = np.sum(L ** 2, axis=0)
        twoE_t = L[0] ** 2 / I1 + L[1] ** 2 / I2 + L[2] ** 2 / I3
        assert np.max(np.abs(M2_t - M2)) / M2 < 1e-13
        assert np.max(np.abs(twoE_t - twoE)) / twoE < 1e-13

    @pytest.mark.parametrize("frac,t_max", [(1e-1, 8000.0), (3e-3, 12000.0)])
    def test_flip_period_across_perturbation_sweep(self, frac, t_max):
        """Closed form tracks the golden across the TC-DZH-03 perturbation range."""
        w0 = np.array([frac * 1e-2, 1e-2, frac * 1e-2])
        M2, twoE = invariants(w0)
        *_, T_LL = ch06_closed_form(M2, twoE)
        ts, w = golden(w0, t_max=t_max)
        meas = flip_interval(ts, w[1])
        assert abs(T_LL / 2.0 - meas) / meas < 1e-8
