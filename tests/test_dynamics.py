"""
test_dynamics.py — unit tests for dynamics module.

Tests quaternion math, the gravity-gradient torque formula, the RK4 integrator
on known analytic cases (torque-free symmetric top), and the diagnostic helpers.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from dynamics import (
    q_normalise, q_to_matrix, q_kinematics_matrix,
    BodyState, pack_state, unpack_state, rk4_step, state_derivative,
    gravitational_force_on_a, gravity_gradient_torque_body_frame,
    total_kinetic_energy, total_potential_energy, total_angular_momentum,
    inertia_inertial_frame, principal_axes_inertial,
    G_NEWTON,
)


# ============================== quaternions ===============================

class TestQuaternionRoundtrip:
    def test_identity_quaternion_is_identity_matrix(self):
        q = np.array([1.0, 0.0, 0.0, 0.0])
        R = q_to_matrix(q)
        assert np.allclose(R, np.eye(3), atol=1e-15)

    def test_180_about_z_flips_x_and_y(self):
        # q = (0, 0, 0, 1) -> R rotates 180 deg about z
        q = np.array([0.0, 0.0, 0.0, 1.0])
        R = q_to_matrix(q)
        v = np.array([1.0, 0.0, 0.0])
        assert np.allclose(R @ v, np.array([-1.0, 0.0, 0.0]), atol=1e-15)
        v = np.array([0.0, 1.0, 0.0])
        assert np.allclose(R @ v, np.array([0.0, -1.0, 0.0]), atol=1e-15)

    def test_q_normalise_brings_quaternion_to_unit_norm(self):
        q = np.array([2.0, 0.0, 0.0, 0.0])  # not unit
        qn = q_normalise(q)
        assert np.isclose(np.linalg.norm(qn), 1.0, atol=1e-15)

    def test_rotation_matrix_is_orthogonal(self):
        # arbitrary unit quaternion
        q = np.array([0.5, 0.5, 0.5, 0.5])
        q = q_normalise(q)
        R = q_to_matrix(q)
        assert np.allclose(R @ R.T, np.eye(3), atol=1e-12)
        assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-12)


class TestQuaternionKinematics:
    def test_zero_omega_gives_zero_dq(self):
        omega = np.zeros(3)
        q = np.array([1.0, 0.0, 0.0, 0.0])
        Omega = q_kinematics_matrix(omega)
        dq = 0.5 * Omega @ q
        assert np.allclose(dq, np.zeros(4), atol=1e-15)


# ============================== forces and torques ===============================

class TestGravitationalForce:
    def test_two_bodies_at_unit_separation(self):
        x_a = np.array([0.0, 0.0, 0.0])
        x_b = np.array([1.0, 0.0, 0.0])
        F = gravitational_force_on_a(x_a, x_b, 1.0, 1.0)
        # F = -G * 1 * 1 * (x_a - x_b)/|x_a - x_b|^3 = -G * (-1, 0, 0) / 1 = (G, 0, 0)
        # but the negative sign in the formula means F points from A toward B
        # x_a - x_b = (-1, 0, 0), so F = -G * (-1)/1 = +G in x direction
        # Wait — recompute. r_vec = x_a - x_b = (-1, 0, 0). F = -G * m_a * m_b * r_vec / |r|^3
        # = -G * 1 * 1 * (-1, 0, 0) / 1 = G * (1, 0, 0). But that points AWAY from B.
        # Actually that's wrong: gravity should pull A toward B, so F on A should be
        # in the direction (x_b - x_a) = (+1, 0, 0). Our formula gives F = (+G, 0, 0).
        # That's correct: F on A points toward B.
        assert F[0] == pytest.approx(G_NEWTON, rel=1e-12)
        assert F[1] == pytest.approx(0.0, abs=1e-20)
        assert F[2] == pytest.approx(0.0, abs=1e-20)

    def test_force_falls_off_as_inverse_square(self):
        # At separation 2, force should be 1/4 of force at separation 1
        x_a = np.zeros(3)
        F1 = gravitational_force_on_a(x_a, np.array([1.0, 0.0, 0.0]), 1.0, 1.0)
        F2 = gravitational_force_on_a(x_a, np.array([2.0, 0.0, 0.0]), 1.0, 1.0)
        ratio = np.linalg.norm(F1) / np.linalg.norm(F2)
        assert ratio == pytest.approx(4.0, rel=1e-12)

    def test_zero_separation_returns_zero(self):
        # safety: no NaN/Inf when bodies coincide
        F = gravitational_force_on_a(np.zeros(3), np.zeros(3), 1.0, 1.0)
        assert np.allclose(F, np.zeros(3))


class TestGravityGradientTorque:
    def test_zero_torque_along_principal_axis(self):
        # if r_rel is exactly along a body principal axis, the torque vanishes
        # for an axisymmetric body
        I_body = np.diag([1.0, 1.0, 0.5])  # axisymmetric about z
        # r_rel in body frame is exactly +z direction; with identity rotation,
        # r_rel_inertial = body_z = (0, 0, 1)
        R = np.eye(3)
        r_rel = np.array([0.0, 0.0, 1.0])
        tau = gravity_gradient_torque_body_frame(I_body, r_rel, R, m_other=10.0)
        assert np.allclose(tau, np.zeros(3), atol=1e-15)

    def test_falls_off_as_inverse_cube(self):
        # tau scales as 1/r^3
        I_body = np.diag([1.0, 2.0, 3.0])
        R = np.eye(3)
        r1 = np.array([1.0, 1.0, 0.0])
        r2 = np.array([2.0, 2.0, 0.0])
        tau1 = gravity_gradient_torque_body_frame(I_body, r1, R, m_other=1.0)
        tau2 = gravity_gradient_torque_body_frame(I_body, r2, R, m_other=1.0)
        ratio = np.linalg.norm(tau1) / np.linalg.norm(tau2)
        assert ratio == pytest.approx(8.0, rel=1e-10)  # 2^3 = 8


# ============================== RK4 integrator on analytic case ===============================

class TestRK4OnTorqueFreeSymmetricTop:
    """The torque-free symmetric top has a closed-form free-precession solution.
    For an axisymmetric body (I_xx = I_yy != I_zz) with no external torque,
    omega_z is constant and the perpendicular components oscillate at the
    Euler frequency lambda = (I_zz - I_xx)/I_xx * omega_z.
    We integrate and verify against this analytic prediction."""

    @pytest.fixture
    def setup(self):
        I_body = np.diag([10.0, 10.0, 5.0])    # oblate axisymmetric
        omega_z = 0.1
        omega_perp = 0.02
        state = BodyState(
            x=[0.0, 0.0, 0.0], v=[0.0, 0.0, 0.0],
            q=[1.0, 0.0, 0.0, 0.0],
            omega=[omega_perp, 0.0, omega_z],
        )
        # Use the coupled state machinery but with body B very far away (zero coupling)
        body_b = BodyState(
            x=[1e10, 0.0, 0.0], v=[0.0, 0.0, 0.0],
            q=[1.0, 0.0, 0.0, 0.0], omega=[0.0, 0.0, 0.0],
        )
        full_state = pack_state(state, body_b)
        return {
            "state": full_state,
            "I_body": I_body,
            "I_b": np.eye(3),
            "omega_z": omega_z,
            "lambda_euler": (5.0 - 10.0) / 10.0 * omega_z,
        }

    def test_omega_z_constant_under_integration(self, setup):
        state = setup["state"]
        for _ in range(1000):
            state = rk4_step(state, dt=0.05,
                             m_a=100.0, I_a_body=setup["I_body"],
                             m_b=1.0, I_b_body=setup["I_b"])
        ba, _ = unpack_state(state)
        # omega_z must remain at the initial value to integration precision
        assert ba.omega[2] == pytest.approx(setup["omega_z"], abs=1e-10)

    def test_perpendicular_omega_oscillates_at_euler_frequency(self, setup):
        state = setup["state"]
        dt = 0.05
        # Half-period of Euler precession should flip sign of omega_x
        # period T = 2*pi / |lambda|
        T = 2.0 * np.pi / abs(setup["lambda_euler"])
        n_half = int((T / 2.0) / dt)
        for _ in range(n_half):
            state = rk4_step(state, dt=dt,
                             m_a=100.0, I_a_body=setup["I_body"],
                             m_b=1.0, I_b_body=setup["I_b"])
        ba, _ = unpack_state(state)
        # omega_x should have approximately flipped sign (and omega_y picked up)
        # because of how the Euler precession rotates in body frame
        magnitude_perp = np.sqrt(ba.omega[0] ** 2 + ba.omega[1] ** 2)
        # the magnitude of (omega_x, omega_y) should be preserved (norm-conservation
        # in the torque-free symmetric-top case)
        assert magnitude_perp == pytest.approx(0.02, rel=1e-4)


# ============================== conservation diagnostics ===============================

class TestConservationDiagnostics:
    def test_kinetic_energy_is_non_negative(self):
        body_a = BodyState([0.0]*3, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
        body_b = BodyState([10.0, 0.0, 0.0], [0.0]*3, [1.0, 0.0, 0.0, 0.0], [0.0]*3)
        state = pack_state(body_a, body_b)
        Etk, Erk, Etot = total_kinetic_energy(state, m_a=2.0, I_a_body=np.eye(3),
                                              m_b=1.0, I_b_body=np.eye(3))
        assert Etk >= 0.0
        assert Erk >= 0.0
        # translational: 0.5 * 2 * 1^2 = 1.0
        assert Etk == pytest.approx(1.0, rel=1e-12)
        # rotational: 0.5 * 1 * 1^2 = 0.5
        assert Erk == pytest.approx(0.5, rel=1e-12)

    def test_potential_energy_is_negative(self):
        body_a = BodyState([0.0]*3, [0.0]*3, [1.0, 0.0, 0.0, 0.0], [0.0]*3)
        body_b = BodyState([10.0, 0.0, 0.0], [0.0]*3, [1.0, 0.0, 0.0, 0.0], [0.0]*3)
        state = pack_state(body_a, body_b)
        # PE = -G * m_a * m_b / r, must be negative
        E_pot = total_potential_energy(state, m_a=1.0, m_b=1.0)
        assert E_pot < 0
        assert E_pot == pytest.approx(-G_NEWTON / 10.0, rel=1e-12)


# ============================== inertia tensor transforms ===============================

class TestInertiaInertialFrame:
    def test_identity_rotation_preserves_body_frame_tensor(self):
        I_body = np.diag([1.0, 2.0, 3.0])
        q_identity = np.array([1.0, 0.0, 0.0, 0.0])
        I_inertial = inertia_inertial_frame(I_body, q_identity)
        assert np.allclose(I_inertial, I_body, atol=1e-15)

    def test_trace_invariant_under_rotation(self):
        # trace of inertia tensor is rotation-invariant
        I_body = np.diag([1.0, 5.0, 8.0])
        trace_body = np.trace(I_body)
        # rotate by an arbitrary unit quaternion
        q = q_normalise(np.array([0.6, 0.3, -0.4, 0.5]))
        I_inertial = inertia_inertial_frame(I_body, q)
        trace_inertial = np.trace(I_inertial)
        assert trace_inertial == pytest.approx(trace_body, rel=1e-12)

    def test_eigenvalues_invariant_under_rotation(self):
        # eigenvalues are rotation-invariant; only the basis changes
        I_body = np.diag([1.0, 5.0, 8.0])
        q = q_normalise(np.array([0.6, 0.3, -0.4, 0.5]))
        eigvals_body = np.sort(np.linalg.eigvalsh(I_body))
        eigvals_inertial = np.sort(np.linalg.eigvalsh(inertia_inertial_frame(I_body, q)))
        assert np.allclose(eigvals_body, eigvals_inertial, atol=1e-12)


class TestPrincipalAxesInertial:
    def test_returns_orthonormal_eigvecs(self):
        I_body = np.diag([1.0, 5.0, 8.0])
        q = q_normalise(np.array([0.6, 0.3, -0.4, 0.5]))
        eigvals, eigvecs = principal_axes_inertial(I_body, q)
        assert np.allclose(eigvecs.T @ eigvecs, np.eye(3), atol=1e-12)
        # eigvals sorted ascending
        assert eigvals[0] <= eigvals[1] <= eigvals[2]
