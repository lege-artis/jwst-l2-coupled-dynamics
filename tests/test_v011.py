"""
test_v011.py — v0.1.1 engineer-tier tests for OQ-4.5 + OQ-3.4 patches.

Two new test classes:
  TestNewtonThirdLawForceBalance  — AG-OQ-4.5-CLOSED
  TestQuadrupoleTidalEnergy       — AG-OQ-3.4-CLOSED
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from dynamics import (
    BodyState, pack_state, state_derivative, total_potential_energy,
    back_reaction_force_one_body, q_to_matrix, q_normalise,
    G_NEWTON,
)


# ============================== Newton's-3rd-law force-balance ===============


class TestNewtonThirdLawForceBalance:
    """AG-OQ-4.5-CLOSED: sum of forces (Kepler + back-reaction) on A + B = 0.

    With back-reaction implemented, m_A * dv_A + m_B * dv_B = 0 to floating-
    point precision because every contributing potential V(rho) depends on
    the relative coordinate only, so CoM momentum is exactly conserved
    (Ch 04 §4.6.4 Newton's-third-law identity).
    """

    @pytest.fixture
    def coupled_state_asymmetric(self):
        """Non-trivial state: both bodies tilted, non-principal-axis spin."""
        m_a = 6173.0
        m_b = 37.0
        I_a = np.diag([23323.0, 23323.0, 15384.0])   # JWST-like prolate
        I_b = np.diag([12.0, 12.0, 70.0])             # probe-like oblate

        # Body A: tilted 0.3 rad about (1,1,0)/sqrt(2)
        theta = 0.3
        ax = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
        q_a = np.concatenate([[np.cos(theta / 2.0)], np.sin(theta / 2.0) * ax])

        ba = BodyState(x=[0.0, 0.0, 0.0], v=[0.0, 0.0, 0.0],
                       q=q_a, omega=[0.02, 0.0, 0.08])
        bb = BodyState(x=[50.0, 0.0, 0.0], v=[0.0, 0.0, 0.0],
                       q=[1.0, 0.0, 0.0, 0.0], omega=[0.0, 0.03, 0.08])
        return pack_state(ba, bb), m_a, I_a, m_b, I_b

    def test_com_force_balance_nominal_separation(self, coupled_state_asymmetric):
        """m_A*dv_A + m_B*dv_B = 0 at 50 m separation (first-cut configuration)."""
        state, m_a, I_a, m_b, I_b = coupled_state_asymmetric
        deriv = state_derivative(state, m_a, I_a, m_b, I_b)
        dv_a = deriv[3:6]
        dv_b = deriv[16:19]
        total_force = m_a * dv_a + m_b * dv_b
        # Should be zero to floating-point precision (internal forces cancel)
        assert np.allclose(total_force, np.zeros(3), atol=1e-20), (
            f"CoM force not zero: {total_force}")

    def test_com_force_balance_close_approach(self):
        """m_A*dv_A + m_B*dv_B = 0 at 10 m separation (back-reaction ~1e-15 larger)."""
        m_a, m_b = 100.0, 10.0
        I_a = np.diag([5.0, 8.0, 12.0])
        I_b = np.diag([1.0, 2.0, 1.5])

        # Non-trivial orientation: 45 deg rotation about z for body A
        q_a = np.array([np.cos(np.pi / 8.0), 0.0, 0.0, np.sin(np.pi / 8.0)])
        q_b = q_normalise(np.array([0.6, 0.3, -0.5, 0.2]))

        ba = BodyState(x=[0.0, 0.0, 0.0], v=[0.0, 0.0, 0.0],
                       q=q_a, omega=[0.1, 0.0, 0.5])
        bb = BodyState(x=[10.0, 0.0, 0.0], v=[0.0, 0.0, 0.0],
                       q=q_b, omega=[0.0, 0.1, 0.3])
        state = pack_state(ba, bb)

        deriv = state_derivative(state, m_a, I_a, m_b, I_b)
        dv_a = deriv[3:6]
        dv_b = deriv[16:19]
        total_force = m_a * dv_a + m_b * dv_b
        assert np.allclose(total_force, np.zeros(3), atol=1e-20), (
            f"CoM force not zero at close approach: {total_force}")

    def test_back_reaction_force_vanishes_for_spherical_body(self):
        """F_back-reaction = 0 for spherical body (§4.6.3 cancellation)."""
        I_sphere = 5.0 * np.eye(3)          # isotropic: I = I_0 * I_3
        r_hat = np.array([1.0, 0.0, 0.0])
        R = np.eye(3)                         # identity rotation
        F_br = back_reaction_force_one_body(I_sphere, r_hat, R,
                                             m_other=100.0, r_mag=50.0)
        assert np.allclose(F_br, np.zeros(3), atol=1e-25), (
            f"Back-reaction not zero for spherical body: {F_br}")

    def test_back_reaction_force_scales_as_r_minus_4(self):
        """Back-reaction force magnitude scales as 1/r^4 (Ch 04 §4.6.2)."""
        I_body = np.diag([3.0, 5.0, 8.0])
        q = q_normalise(np.array([0.7, 0.1, -0.3, 0.4]))
        R = q_to_matrix(q)
        r_hat = np.array([1.0, 0.0, 0.0])

        F1 = back_reaction_force_one_body(I_body, r_hat, R, m_other=50.0, r_mag=1.0)
        F2 = back_reaction_force_one_body(I_body, r_hat, R, m_other=50.0, r_mag=2.0)
        ratio = np.linalg.norm(F1) / np.linalg.norm(F2)
        # doubling r should reduce force by factor 2^4 = 16
        assert ratio == pytest.approx(16.0, rel=1e-10)


# ============================== Quadrupole-tidal energy ======================


class TestQuadrupoleTidalEnergy:
    """AG-OQ-3.4-CLOSED: total_potential_energy includes quadrupole-tidal term.

    Verifies the Ch 03 §3.4.3 boxed augmentation:
      V_tidal_X = -(G m_other / (2 r^3)) * (tr(I_X) - 3 r_hat_body . I_X . r_hat_body)
    """

    @pytest.fixture
    def aligned_state(self):
        """Bodies aligned along x-axis with identity rotation matrices."""
        m_a, m_b = 100.0, 50.0
        I_a = np.diag([2.0, 3.0, 5.0])
        I_b = np.diag([1.0, 1.5, 4.0])
        ba = BodyState(x=[0.0, 0.0, 0.0], v=[0.0, 0.0, 0.0],
                       q=[1.0, 0.0, 0.0, 0.0], omega=[0.0, 0.0, 0.0])
        bb = BodyState(x=[100.0, 0.0, 0.0], v=[0.0, 0.0, 0.0],
                       q=[1.0, 0.0, 0.0, 0.0], omega=[0.0, 0.0, 0.0])
        return pack_state(ba, bb), m_a, m_b, I_a, I_b

    def test_tidal_energy_scales_linearly_with_inertia_scale(self, aligned_state):
        """V_tidal scales linearly with I: scaling I by alpha scales V_tidal by alpha."""
        state, m_a, m_b, I_a, I_b = aligned_state
        alpha = 3.7

        V1 = total_potential_energy(state, m_a, m_b, I_a, I_b)
        V_kepler = total_potential_energy(state, m_a, m_b)   # no inertia args
        dV1 = V1 - V_kepler                                   # tidal part only

        V2 = total_potential_energy(state, m_a, m_b, alpha * I_a, alpha * I_b)
        dV2 = V2 - V_kepler

        # tidal contribution must scale with alpha
        assert dV2 == pytest.approx(alpha * dV1, rel=1e-12)

    def test_tidal_energy_zero_for_spherical_bodies(self, aligned_state):
        """V_tidal = 0 for spherical inertia tensors (§3.8.2 decoupling)."""
        state, m_a, m_b, _, _ = aligned_state
        I_sphere_a = 4.0 * np.eye(3)
        I_sphere_b = 2.5 * np.eye(3)
        V_with = total_potential_energy(state, m_a, m_b, I_sphere_a, I_sphere_b)
        V_kepler = total_potential_energy(state, m_a, m_b)
        assert V_with == pytest.approx(V_kepler, rel=1e-14)

    def test_tidal_energy_correct_analytic_value_aligned_axis(self, aligned_state):
        """Manually verify V_tidal_A for body A aligned with separation axis (+x).

        With R = I (identity), r_hat_body_A = (1, 0, 0), so:
          Phi_A = tr(I_A) - 3 * I_A[0,0]
          V_tidal_A = -(G * m_B / (2 * r^3)) * Phi_A
        """
        state, m_a, m_b, I_a, I_b = aligned_state
        r = 100.0
        # body A: I_a = diag(2, 3, 5), r_hat_body = (1, 0, 0)
        Phi_A_expected = np.trace(I_a) - 3.0 * I_a[0, 0]   # 10 - 6 = 4
        V_tidal_A_expected = -(G_NEWTON * m_b / (2.0 * r**3)) * Phi_A_expected
        # body B: I_b = diag(1, 1.5, 4), r_hat_body_B = R_B^T @ (-r_hat)
        # r_hat = (1,0,0), R_B = I, so r_hat_body_B = (-1,0,0)
        # Phi_B = tr(I_b) - 3 * ((-1)^2 * I_b[0,0]) = tr(I_b) - 3 * I_b[0,0]
        Phi_B_expected = np.trace(I_b) - 3.0 * I_b[0, 0]   # 6.5 - 3 = 3.5
        V_tidal_B_expected = -(G_NEWTON * m_a / (2.0 * r**3)) * Phi_B_expected

        V_total_expected = -(G_NEWTON * m_a * m_b / r) + V_tidal_A_expected + V_tidal_B_expected
        V_computed = total_potential_energy(state, m_a, m_b, I_a, I_b)
        assert V_computed == pytest.approx(V_total_expected, rel=1e-12)

    def test_tidal_energy_independent_of_separation_sign(self):
        """V_tidal is even in r_hat (Ch 03 §3.4.3 — quadratic in direction)."""
        m_a, m_b = 50.0, 30.0
        I_a = np.diag([2.0, 4.0, 6.0])
        I_b = np.diag([1.0, 2.0, 3.0])
        # arbitrary rotation for body A
        q_a = q_normalise(np.array([0.5, 0.3, -0.6, 0.4]))

        ba1 = BodyState(x=[0.0, 0.0, 0.0], v=[0.0]*3, q=q_a, omega=[0.0]*3)
        bb1 = BodyState(x=[0.0, 80.0, 0.0], v=[0.0]*3,
                        q=[1.0, 0.0, 0.0, 0.0], omega=[0.0]*3)
        ba2 = BodyState(x=[0.0, 80.0, 0.0], v=[0.0]*3, q=q_a, omega=[0.0]*3)
        bb2 = BodyState(x=[0.0, 0.0, 0.0], v=[0.0]*3,
                        q=[1.0, 0.0, 0.0, 0.0], omega=[0.0]*3)

        V1 = total_potential_energy(pack_state(ba1, bb1), m_a, m_b, I_a, I_b)
        V2 = total_potential_energy(pack_state(ba2, bb2), m_b, m_a, I_b, I_a)
        # By symmetry (swap A and B), both should give the same total V
        assert V1 == pytest.approx(V2, rel=1e-12)

    def test_tidal_energy_reduces_to_kepler_when_no_inertia_passed(self, aligned_state):
        """Backward-compatible: omitting inertia args returns Kepler-only PE."""
        state, m_a, m_b, I_a, I_b = aligned_state
        V_compat = total_potential_energy(state, m_a, m_b)
        V_kepler = -G_NEWTON * m_a * m_b / 100.0
        assert V_compat == pytest.approx(V_kepler, rel=1e-14)
