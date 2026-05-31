"""
dynamics.py — coupled translational + rotational dynamics for two rigid bodies
under mutual gravitational attraction with gravity-gradient torques.

State per body: x (3), v (3), q (4 unit-quaternion), omega (3 in body frame)
Total state vector for two bodies: 13 + 13 = 26 components.

Conventions:
- Quaternion convention: q = (qw, qx, qy, qz) with qw = scalar part
- Rotation matrix R = R(q) maps body-frame vectors to inertial-frame vectors
- Angular velocity omega is in the BODY frame (Euler's equations canonical form)
- Inertia tensor I_body is constant in the body frame; transform to inertial via R . I . R^T

Equations of motion (v0.1.1: includes back-reaction force + tidal potential):
  Translation (Newton, inertial frame):
    dx/dt = v
    dv/dt = (F_kepler + F_back-reaction) / m
    F_kepler = -G * m_self * m_other * (x_self - x_other) / |x_self - x_other|^3
    F_back-reaction: Ch 04 §4.6.2 boxed (OQ-4.5)

  Rotation (Euler, body frame):
    I_body . d(omega)/dt + omega x (I_body . omega) = tau_body
    dq/dt = 0.5 * Omega(omega) . q  (quaternion kinematics)

Gravity-gradient torque on a body at separation r from a point mass m_other:
  tau_body = (3 * G * m_other / |r|^3) * (r_hat_body) x (I_body . r_hat_body)
  where r_hat_body = R^T . r_hat_inertial (relative direction in body frame)

References:
- Goldstein, Poole, Safko (2002) "Classical Mechanics" 3rd ed §5.6 (Euler eq)
- Hughes (2004) "Spacecraft Attitude Dynamics" §3, §5 (gravity-gradient torque)
- Markley & Crassidis (2014) "Fundamentals of Spacecraft Attitude Determination" §3
"""
from __future__ import annotations
import numpy as np

# --- physical constants ---
G_NEWTON = 6.67430e-11  # m^3 / (kg s^2) — gravitational constant (CODATA 2018)


# --- quaternion utilities ---

def q_normalise(q: np.ndarray) -> np.ndarray:
    """Renormalise a quaternion to unit length (counters integration drift)."""
    return q / np.linalg.norm(q)


def q_to_matrix(q: np.ndarray) -> np.ndarray:
    """Rotation matrix from unit quaternion (qw, qx, qy, qz).
       R maps body-frame vectors to inertial-frame vectors: v_inertial = R . v_body."""
    qw, qx, qy, qz = q
    R = np.array([
        [1.0 - 2.0*(qy*qy + qz*qz), 2.0*(qx*qy - qz*qw),     2.0*(qx*qz + qy*qw)],
        [2.0*(qx*qy + qz*qw),       1.0 - 2.0*(qx*qx + qz*qz), 2.0*(qy*qz - qx*qw)],
        [2.0*(qx*qz - qy*qw),       2.0*(qy*qz + qx*qw),     1.0 - 2.0*(qx*qx + qy*qy)],
    ])
    return R


def q_kinematics_matrix(omega: np.ndarray) -> np.ndarray:
    """Quaternion rate matrix Omega(omega) such that dq/dt = 0.5 * Omega(omega) . q.
       omega is in the BODY frame."""
    wx, wy, wz = omega
    return np.array([
        [ 0.0, -wx, -wy, -wz],
        [ wx,  0.0,  wz, -wy],
        [ wy, -wz,  0.0,  wx],
        [ wz,  wy, -wx,  0.0],
    ])


# --- body-state dataclass ---

class BodyState:
    """State of one rigid body: position, velocity, attitude, angular velocity."""
    __slots__ = ("x", "v", "q", "omega")

    def __init__(self, x: np.ndarray, v: np.ndarray,
                 q: np.ndarray, omega: np.ndarray):
        self.x = np.array(x, dtype=float)
        self.v = np.array(v, dtype=float)
        self.q = np.array(q, dtype=float)
        self.omega = np.array(omega, dtype=float)

    def as_flat(self) -> np.ndarray:
        """Flatten the 13-component state to a 1D array."""
        return np.concatenate([self.x, self.v, self.q, self.omega])

    @classmethod
    def from_flat(cls, arr: np.ndarray) -> "BodyState":
        return cls(arr[0:3], arr[3:6], arr[6:10], arr[10:13])


def pack_state(body_a: BodyState, body_b: BodyState) -> np.ndarray:
    return np.concatenate([body_a.as_flat(), body_b.as_flat()])


def unpack_state(state: np.ndarray) -> tuple[BodyState, BodyState]:
    return BodyState.from_flat(state[:13]), BodyState.from_flat(state[13:26])


# --- force + torque ---

def gravitational_force_on_a(x_a: np.ndarray, x_b: np.ndarray,
                              m_a: float, m_b: float) -> np.ndarray:
    """Newtonian gravitational force on body A due to body B (inertial frame)."""
    r_vec = x_a - x_b
    r_mag = np.linalg.norm(r_vec)
    if r_mag < 1e-12:
        return np.zeros(3)
    return -G_NEWTON * m_a * m_b * r_vec / (r_mag ** 3)


def gravity_gradient_torque_body_frame(I_body: np.ndarray,
                                         r_rel_inertial: np.ndarray,
                                         R_body_to_inertial: np.ndarray,
                                         m_other: float) -> np.ndarray:
    """Gravity-gradient torque on a body of inertia I_body due to a point mass
       m_other at relative inertial position r_rel = (x_other - x_self).

       Returned vector is in the BODY frame (Euler's equations canonical form).

       Formula (Hughes 2004 §5.3): tau_body = (3 G m_other / r^3) * r_hat_body x (I . r_hat_body)
       where r_hat_body = R^T . r_hat_inertial points from self to other in body frame."""
    r_mag = np.linalg.norm(r_rel_inertial)
    if r_mag < 1e-12:
        return np.zeros(3)
    r_hat_inertial = r_rel_inertial / r_mag
    r_hat_body = R_body_to_inertial.T @ r_hat_inertial
    I_rhat = I_body @ r_hat_body
    tau = (3.0 * G_NEWTON * m_other / (r_mag ** 3)) * np.cross(r_hat_body, I_rhat)
    return tau


def back_reaction_force_one_body(I_body: np.ndarray,
                                  r_hat: np.ndarray,
                                  R: np.ndarray,
                                  m_other: float,
                                  r_mag: float) -> np.ndarray:
    """Back-reaction force on the relative orbit from one body's tidal coupling.

    Ch 04 §4.6.2 boxed (OQ-4.5):
      F = -(3Gm/(2r^4)) Phi r_hat
          - (3Gm/r^4) (I3 - r_hat r_hat^T) R (I_body . r_hat_body)
    where r_hat_body = R^T . r_hat  and  Phi = tr(I) - 3 r_hat_body . I . r_hat_body.

    r_hat is always the unit vector from body A to body B in the inertial frame,
    used identically for both body A and body B contributions (§4.6.2, §4.6.4).
    The total back-reaction force F_br_A + F_br_B enters the relative equation
    mu*rho_ddot = F_kepler + F_br_total.  Force on body A = -F_br_total;
    force on body B = +F_br_total (Newton's third law, §4.6.4)."""
    r_hat_body = R.T @ r_hat
    Phi = np.trace(I_body) - 3.0 * (r_hat_body @ I_body @ r_hat_body)
    coeff = 3.0 * G_NEWTON * m_other / r_mag**4
    I_rhat_inertial = R @ (I_body @ r_hat_body)
    # (I3 - r_hat r_hat^T) R (I_body r_hat_body) = transverse projection
    transverse = I_rhat_inertial - r_hat * (r_hat @ I_rhat_inertial)
    return -0.5 * coeff * Phi * r_hat - coeff * transverse


# --- derivatives: dstate/dt ---

def state_derivative(state: np.ndarray,
                      m_a: float, I_a_body: np.ndarray,
                      m_b: float, I_b_body: np.ndarray) -> np.ndarray:
    """Compute dstate/dt for the 26-component coupled-body state."""
    body_a, body_b = unpack_state(state)

    # Separation vector and rotation matrices (used by both translation and rotation)
    r_vec = body_b.x - body_a.x          # vector from A to B in inertial frame
    r_mag = np.linalg.norm(r_vec)
    R_a = q_to_matrix(body_a.q)
    R_b = q_to_matrix(body_b.q)

    # Kepler translational dynamics
    F_a = gravitational_force_on_a(body_a.x, body_b.x, m_a, m_b)
    F_b = -F_a   # Newton's 3rd law (Kepler piece)
    dx_a = body_a.v
    dv_a = F_a / m_a
    dx_b = body_b.v
    dv_b = F_b / m_b

    # Back-reaction force (Ch 04 §4.6.2, OQ-4.5 CLOSED)
    # Both bodies use the same r_hat = (x_B - x_A)/r (§4.6.2 symmetry).
    # Force on A = -(F_br_A + F_br_B); force on B = +(F_br_A + F_br_B).
    if r_mag >= 1e-12:
        r_hat = r_vec / r_mag
        F_br_a = back_reaction_force_one_body(I_a_body, r_hat, R_a, m_b, r_mag)
        F_br_b = back_reaction_force_one_body(I_b_body, r_hat, R_b, m_a, r_mag)
        F_br_total = F_br_a + F_br_b
        dv_a -= F_br_total / m_a
        dv_b += F_br_total / m_b

    # Rotational dynamics — gravity-gradient torque from the other body
    tau_a = gravity_gradient_torque_body_frame(I_a_body, r_vec, R_a, m_b)
    tau_b = gravity_gradient_torque_body_frame(I_b_body, -r_vec, R_b, m_a)

    # Euler's equations: I . d(omega)/dt = tau - omega x (I . omega)
    I_omega_a = I_a_body @ body_a.omega
    I_omega_b = I_b_body @ body_b.omega
    domega_a = np.linalg.solve(I_a_body, tau_a - np.cross(body_a.omega, I_omega_a))
    domega_b = np.linalg.solve(I_b_body, tau_b - np.cross(body_b.omega, I_omega_b))

    # Quaternion kinematics: dq/dt = 0.5 * Omega(omega) . q
    dq_a = 0.5 * q_kinematics_matrix(body_a.omega) @ body_a.q
    dq_b = 0.5 * q_kinematics_matrix(body_b.omega) @ body_b.q

    return np.concatenate([dx_a, dv_a, dq_a, domega_a,
                            dx_b, dv_b, dq_b, domega_b])


# --- RK4 integrator ---

def rk4_step(state: np.ndarray, dt: float,
             m_a: float, I_a_body: np.ndarray,
             m_b: float, I_b_body: np.ndarray) -> np.ndarray:
    """One step of the classical 4th-order Runge-Kutta integrator.
       Reference: Press et al. (2007) Numerical Recipes 3rd ed §17.1, eq. 17.1.3.
       After the step, the quaternions are renormalised to combat integration drift."""
    k1 = state_derivative(state, m_a, I_a_body, m_b, I_b_body)
    k2 = state_derivative(state + 0.5 * dt * k1, m_a, I_a_body, m_b, I_b_body)
    k3 = state_derivative(state + 0.5 * dt * k2, m_a, I_a_body, m_b, I_b_body)
    k4 = state_derivative(state + dt * k3, m_a, I_a_body, m_b, I_b_body)
    new_state = state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    # Renormalise quaternions in-place
    new_state[6:10] = q_normalise(new_state[6:10])
    new_state[19:23] = q_normalise(new_state[19:23])
    return new_state


# --- diagnostics ---

def total_kinetic_energy(state: np.ndarray,
                          m_a: float, I_a_body: np.ndarray,
                          m_b: float, I_b_body: np.ndarray) -> tuple[float, float, float]:
    """Returns (E_trans_total, E_rot_total, E_trans + E_rot)."""
    body_a, body_b = unpack_state(state)
    E_trans = 0.5 * m_a * np.dot(body_a.v, body_a.v) + 0.5 * m_b * np.dot(body_b.v, body_b.v)
    E_rot = 0.5 * body_a.omega @ I_a_body @ body_a.omega
    E_rot += 0.5 * body_b.omega @ I_b_body @ body_b.omega
    return E_trans, E_rot, E_trans + E_rot


def total_potential_energy(state: np.ndarray, m_a: float, m_b: float,
                            I_a_body: np.ndarray | None = None,
                            I_b_body: np.ndarray | None = None) -> float:
    """Gravitational PE: monopole-monopole (Kepler) + optional quadrupole-tidal terms.

    When I_a_body and I_b_body are supplied, adds the Ch 03 §3.4.3 tidal pieces
    V_tidal_A + V_tidal_B so the conservation diagnostic tracks the full
    quadrupole-truncated V rather than the monopole-only term (OQ-3.4 CLOSED).
    Omitting the inertia tensors returns the v0.1.0 monopole-only result."""
    body_a, body_b = unpack_state(state)
    r_vec = body_b.x - body_a.x          # from A to B
    r_mag = np.linalg.norm(r_vec)
    if r_mag < 1e-12:
        return 0.0
    V_kepler = -G_NEWTON * m_a * m_b / r_mag
    if I_a_body is None or I_b_body is None:
        return V_kepler
    # Quadrupole-tidal augmentation (Ch 03 §3.4.3 boxed, OQ-3.4)
    r_hat = r_vec / r_mag
    R_a = q_to_matrix(body_a.q)
    R_b = q_to_matrix(body_b.q)
    r_hat_body_a = R_a.T @ r_hat
    V_tidal_a = -(G_NEWTON * m_b / (2.0 * r_mag**3)) * (
        np.trace(I_a_body) - 3.0 * (r_hat_body_a @ I_a_body @ r_hat_body_a))
    r_hat_body_b = R_b.T @ (-r_hat)      # direction from B toward A in B's body frame
    V_tidal_b = -(G_NEWTON * m_a / (2.0 * r_mag**3)) * (
        np.trace(I_b_body) - 3.0 * (r_hat_body_b @ I_b_body @ r_hat_body_b))
    return V_kepler + V_tidal_a + V_tidal_b


def total_angular_momentum(state: np.ndarray,
                            m_a: float, I_a_body: np.ndarray,
                            m_b: float, I_b_body: np.ndarray) -> np.ndarray:
    """Total system angular momentum in the inertial frame.
       L_total = L_orbital_A + L_orbital_B + L_spin_A + L_spin_B
       where L_spin_X = R_X . (I_X_body . omega_X)."""
    body_a, body_b = unpack_state(state)

    L_orbital_a = m_a * np.cross(body_a.x, body_a.v)
    L_orbital_b = m_b * np.cross(body_b.x, body_b.v)

    R_a = q_to_matrix(body_a.q)
    R_b = q_to_matrix(body_b.q)
    L_spin_a = R_a @ (I_a_body @ body_a.omega)
    L_spin_b = R_b @ (I_b_body @ body_b.omega)

    return L_orbital_a + L_orbital_b + L_spin_a + L_spin_b


def inertia_inertial_frame(I_body: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Transform inertia tensor from body frame to inertial frame:
       I_inertial = R . I_body . R^T."""
    R = q_to_matrix(q)
    return R @ I_body @ R.T


def principal_axes_inertial(I_body: np.ndarray, q: np.ndarray
                             ) -> tuple[np.ndarray, np.ndarray]:
    """Eigendecomposition of the inertia tensor expressed in the inertial frame.
       Returns (eigenvalues_ascending, eigvecs_inertial[:, k]) where each
       eigenvector is a principal axis direction in the inertial frame."""
    I_inertial = inertia_inertial_frame(I_body, q)
    return np.linalg.eigh(I_inertial)
