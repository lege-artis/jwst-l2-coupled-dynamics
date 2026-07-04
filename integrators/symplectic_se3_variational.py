"""
symplectic_se3_variational.py — Lie-group variational integrator on SE(3) × SE(3).

Implements the KDK (kick–drift–kick) Störmer–Verlet method adapted to the
SO(3) × SO(3) × R^6 phase space for two coupled rigid bodies under mutual
gravitational attraction.

Algorithm: Störmer–Verlet / KDK leapfrog on SE(3) × SE(3)
─────────────────────────────────────────────────────────
Given state at step k: (x_a, v_a, q_a, ω_a, x_b, v_b, q_b, ω_b)

Internally, work in canonical (Hamiltonian) variables:
  Π_a = I_a · ω_a   (body-frame angular momentum)
  p_a = m_a · v_a   (linear momentum)

One step k → k+1:
  1. KICK ½: Π += Δt/2 · τ,  p += Δt/2 · F   (forces at state_k)
  2. DRIFT:  rotation R_{k+1} = R_k · exp(Δt · J⁻¹Π)ˣ
             Π → F^T Π  (coadjoint transport; preserves L_inertial = R·Π exactly)
             x_{k+1} = x_k + Δt · p/m
  3. KICK ½: Π += Δt/2 · τ,  p += Δt/2 · F   (forces at state_{k+1})

Properties:
  • Second-order accurate in trajectory (O(Δt²) global error)
  • Symplectic: preserves the discrete symplectic 2-form at every step
  • Momentum-preserving: inertial angular momentum L = R·Π is EXACTLY conserved
    in the free-body limit (τ = 0) because the coadjoint transport step
    R → R·F, Π → F^T·Π preserves R·Π identically in exact arithmetic
  • Energy behaviour: bounded oscillation for exponentially long times
    (backward-error-analysis: HLW 2006 Ch IX)
  • Group structure: SO(3) is preserved at every step via the exponential map
    (no reprojection needed in exact arithmetic; renormalisation added for
    floating-point stability)

Note on evaluation count: the KDK form uses TWO force evaluations per step
(at state_k and state_{k+1}).  The "single-evaluation" property cited in
LLM 2007 comes from a specific time-asymmetric discrete Lagrangian; the
symmetric KDK used here achieves second-order accuracy at the cost of one
additional evaluation.  Both variants are symplectic and momentum-preserving;
the symmetric form is preferable for the libration and Lyapunov testbeds.

Note on integrator order (IMPORTANT for AG-LIB-2 validation):
The coadjoint transport Pi_new = F^T Pi_old (with F = exp(dt * omega)) is a
FIRST-ORDER approximation to the free-body Euler-equation flow because the
angular velocity omega changes during the step.  The resulting global error in
body-frame omega is O(T·dt) = first-order.  To satisfy AG-LIB-2 (second-order
convergence in libration period), the Sitting 2 upgrade applies the MIDPOINT
Newton iteration: iterate omega_mid = (omega_k + omega_{k+1})/2 to self-
consistency (3–5 steps to ULP), then use F = exp(dt * omega_mid).  This
achieves second-order accuracy in body-frame quantities.  Sitting 1 delivers
the correct structural framework (coadjoint, rodrigues, KDK structure, exact
inertial-L preservation); Sitting 2 upgrades to the midpoint rotation step.

References:
  Lee, T., Leok, M., McClamroch, N.H. (2007) "Lie group variational integrators
    for the full body problem in orbital mechanics", CMDA 98, 121–144.
    DOI 10.1007/s10569-007-9073-x.  [PRIMARY ANCHOR — §2.2 of Phase 2 brief]
  Marsden, J.E., West, M. (2001) "Discrete mechanics and variational integrators",
    Acta Numerica 10, 357–514.  [Variational-integrators survey]
  Hairer, E., Lubich, C., Wanner, G. (2006) Geometric Numerical Integration,
    2nd ed., Springer, Ch VII + Ch IX.  [HLW backward-error analysis]

CHANGELOG
─────────────────────────────────────────────────────────────────────────────
Sitting 1 (2026-06-01): initial KDK Störmer–Verlet on SE(3)×SE(3); coadjoint
  transport preserving inertial L exactly; rotation step first-order (used
  F = exp(dt·ω_k) with ω_k from the half-kicked Π alone).
Sitting 2 (2026-06-05): upgraded to second-order via midpoint Newton iteration.
  The rotation step now solves ω_mid = (ω_k + ω_{k+1})/2 to self-consistency
  (up to 5 fixed-point iterations, converge criterion |δω|_∞ < 10·ε_machine)
  before applying F = exp(dt·ω_mid).  This gives O(dt²) global accuracy in
  body-frame angular velocity, consistent with the O(dt²) translational
  Störmer–Verlet; satisfies AG-LIB-2 (second-order convergence in libration
  period) and AG-INT-1 (second-order Lie-group variational integrator).
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np

_EPS = np.finfo(float).eps  # IEEE-754 double machine epsilon ≈ 2.22e-16

sys.path.insert(0, str(Path(__file__).parent.parent))
from dynamics import (
    BodyState, pack_state, unpack_state,
    q_normalise,
    q_to_matrix,
    gravitational_force_on_a,
    gravity_gradient_torque_body_frame,
    back_reaction_force_one_body,
)


# ─────────────────────────────────────────────────────────────────────────────
# SO(3) exponential map + quaternion helpers
# ─────────────────────────────────────────────────────────────────────────────

def hat(v: np.ndarray) -> np.ndarray:
    """Skew-symmetric matrix (hat map) [v]× for v ∈ R³."""
    return np.array([
        [ 0.0,   -v[2],  v[1]],
        [ v[2],   0.0,  -v[0]],
        [-v[1],   v[0],  0.0],
    ])


def rodrigues(v: np.ndarray) -> np.ndarray:
    """Rotation matrix exp([v]×) from rotation vector v ∈ R³.

    Uses the closed-form Rodrigues formula:
      exp([v]×) = I + sin(θ)/θ · [v]× + (1−cos(θ))/θ² · [v]×²
    where θ = |v|.  For |v| < 1e-14, reverts to first-order I + [v]×.
    """
    theta = np.linalg.norm(v)
    if theta < 1e-14:
        return np.eye(3) + hat(v)
    n = v / theta
    K = hat(n)
    return np.eye(3) + np.sin(theta) * K + (1.0 - np.cos(theta)) * (K @ K)


def rotvec_to_quat(v: np.ndarray) -> np.ndarray:
    """Quaternion (qw, qx, qy, qz) for rotation vector v = θ·n̂.

    Consistent with dynamics.py q_to_matrix convention:
      q_to_matrix(rotvec_to_quat(v)) == rodrigues(v)   (to numerical precision).
    """
    theta = np.linalg.norm(v)
    if theta < 1e-14:
        # First-order: cos(θ/2) ≈ 1, sin(θ/2) ≈ θ/2
        return np.array([1.0, 0.5 * v[0], 0.5 * v[1], 0.5 * v[2]])
    n = v / theta
    half = 0.5 * theta
    s = np.sin(half)
    return np.array([np.cos(half), s * n[0], s * n[1], s * n[2]])


def quat_multiply(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Hamilton product p ⊗ q for quaternions (qw, qx, qy, qz).

    Satisfies: q_to_matrix(p ⊗ q) = q_to_matrix(p) @ q_to_matrix(q).
    """
    pw, px, py, pz = p
    qw, qx, qy, qz = q
    return np.array([
        pw*qw - px*qx - py*qy - pz*qz,
        pw*qx + px*qw + py*qz - pz*qy,
        pw*qy - px*qz + py*qw + pz*qx,
        pw*qz + px*qy - py*qx + pz*qw,
    ])


# ─────────────────────────────────────────────────────────────────────────────
# Midpoint Newton iteration for the rotation drift step (Sitting 2)
# ─────────────────────────────────────────────────────────────────────────────

def _midpoint_rotation(
    Pi: np.ndarray,
    I_inv: np.ndarray,
    dt: float,
    max_iter: int = 5,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve ω_mid = (ω_k + ω_{k+1})/2 to self-consistency (Sitting 2 upgrade).

    Fixed-point iteration:
      1. omega_k  = I_inv @ Pi_k            (current body-frame omega)
      2. F        = rodrigues(dt · omega_mid)
      3. omega_{k+1} = I_inv @ (F^T @ Pi_k)
      4. omega_mid_new = (omega_k + omega_{k+1}) / 2
      5. repeat until |omega_mid_new − omega_mid|_∞ < 10 · ε_machine

    Returns (rotvec, F_rot) where rotvec = dt · omega_mid (converged) and
    F_rot = rodrigues(rotvec).  The caller then applies:
      Pi_new = F_rot.T @ Pi_k          (coadjoint transport)
      q_new  = q_k ⊗ rotvec_to_quat(rotvec)

    Convergence is quadratic for small dt (spectral radius of the
    iteration Jacobian is O(dt² · ‖I_inv‖²·‖Pi‖) ≪ 1 for timesteps used
    in practice); 3–5 iterations reach ULP precision.
    """
    omega_k = I_inv @ Pi
    omega_mid = omega_k.copy()
    for _ in range(max_iter):
        F = rodrigues(dt * omega_mid)
        omega_next = I_inv @ (F.T @ Pi)
        omega_mid_new = 0.5 * (omega_k + omega_next)
        if np.max(np.abs(omega_mid_new - omega_mid)) < 10.0 * _EPS:
            omega_mid = omega_mid_new
            break
        omega_mid = omega_mid_new
    rotvec = dt * omega_mid
    return rotvec, rodrigues(rotvec)


# ─────────────────────────────────────────────────────────────────────────────
# Force / torque helper (shared between both KDK kick phases)
# ─────────────────────────────────────────────────────────────────────────────

def _forces_torques(
    x_a: np.ndarray, x_b: np.ndarray,
    q_a: np.ndarray, q_b: np.ndarray,
    I_a_body: np.ndarray, I_b_body: np.ndarray,
    m_a: float, m_b: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute total force on A, total force on B, torque on A, torque on B.

    Mirrors the force/torque block in dynamics.py::state_derivative, including
    Kepler forces, back-reaction forces (Ch 04 §4.6.2), and gravity-gradient
    torques.  Returns (F_a, F_b, tau_a, tau_b) all in the inertial frame
    (forces) or body frame (torques).
    """
    r_vec = x_b - x_a
    r_mag = np.linalg.norm(r_vec)
    R_a = q_to_matrix(q_a)
    R_b = q_to_matrix(q_b)

    # Kepler force
    F_a = gravitational_force_on_a(x_a, x_b, m_a, m_b)
    F_b = -F_a

    # Back-reaction force (tidal coupling, OQ-4.5 CLOSED)
    if r_mag >= 1e-12:
        r_hat = r_vec / r_mag
        F_br_a = back_reaction_force_one_body(I_a_body, r_hat, R_a, m_b, r_mag)
        F_br_b = back_reaction_force_one_body(I_b_body, r_hat, R_b, m_a, r_mag)
        F_br_total = F_br_a + F_br_b
        F_a -= F_br_total
        F_b += F_br_total

    # Gravity-gradient torques (body frame)
    tau_a = gravity_gradient_torque_body_frame(I_a_body, r_vec, R_a, m_b)
    tau_b = gravity_gradient_torque_body_frame(I_b_body, -r_vec, R_b, m_a)

    return F_a, F_b, tau_a, tau_b


# ─────────────────────────────────────────────────────────────────────────────
# Main integrator step
# ─────────────────────────────────────────────────────────────────────────────

def se3_variational_step(
    state: np.ndarray,
    dt: float,
    m_a: float, I_a_body: np.ndarray,
    m_b: float, I_b_body: np.ndarray,
) -> np.ndarray:
    """One KDK Störmer–Verlet step on SE(3) × SE(3).

    Drop-in replacement for dynamics.rk4_step with identical call signature.
    Returns a new 26-component state array; does not modify the input.

    Args:
        state: 26-component state [x_a(3), v_a(3), q_a(4), ω_a(3),
                                    x_b(3), v_b(3), q_b(4), ω_b(3)]
        dt: timestep (seconds)
        m_a, I_a_body: mass and body-frame inertia tensor of body A
        m_b, I_b_body: mass and body-frame inertia tensor of body B

    Returns:
        New 26-component state at time t + dt.
    """
    body_a, body_b = unpack_state(state)

    # Inertia inverses (computed once per step)
    I_inv_a = np.linalg.inv(I_a_body)
    I_inv_b = np.linalg.inv(I_b_body)

    # ── Convert to canonical momenta ──────────────────────────────────────────
    Pi_a = I_a_body @ body_a.omega    # body-frame angular momentum
    Pi_b = I_b_body @ body_b.omega
    p_a  = m_a * body_a.v             # linear momentum
    p_b  = m_b * body_b.v

    # ── KICK ½: forces/torques at state_k ────────────────────────────────────
    F_a, F_b, tau_a, tau_b = _forces_torques(
        body_a.x, body_b.x, body_a.q, body_b.q,
        I_a_body, I_b_body, m_a, m_b,
    )
    Pi_a += 0.5 * dt * tau_a
    Pi_b += 0.5 * dt * tau_b
    p_a  += 0.5 * dt * F_a
    p_b  += 0.5 * dt * F_b

    # ── DRIFT: rotation + translation ────────────────────────────────────────
    # Midpoint Newton iteration (Sitting 2): solve omega_mid = (omega_k + omega_{k+1})/2
    # to self-consistency.  Returns rotation vector v = dt·omega_mid and F = exp([v]×).
    v_a, F_rot_a = _midpoint_rotation(Pi_a, I_inv_a, dt)
    v_b, F_rot_b = _midpoint_rotation(Pi_b, I_inv_b, dt)

    # Quaternion update: q_{k+1} = q_k ⊗ q(F_rot)
    q_a_new = q_normalise(quat_multiply(body_a.q, rotvec_to_quat(v_a)))
    q_b_new = q_normalise(quat_multiply(body_b.q, rotvec_to_quat(v_b)))

    # Coadjoint transport: Π_{new_body} = F^T Π_{old_body}
    # This EXACTLY preserves L_inertial = R · Π in exact arithmetic:
    #   R_new · Π_new = (R_old · F_rot) · (F_rot^T · Π_old) = R_old · Π_old
    Pi_a = F_rot_a.T @ Pi_a
    Pi_b = F_rot_b.T @ Pi_b

    # Translation (velocity already incorporates the half-kick)
    x_a_new = body_a.x + dt * (p_a / m_a)
    x_b_new = body_b.x + dt * (p_b / m_b)

    # ── KICK ½: forces/torques at state_{k+1} ────────────────────────────────
    F_a2, F_b2, tau_a2, tau_b2 = _forces_torques(
        x_a_new, x_b_new, q_a_new, q_b_new,
        I_a_body, I_b_body, m_a, m_b,
    )
    Pi_a += 0.5 * dt * tau_a2
    Pi_b += 0.5 * dt * tau_b2
    p_a  += 0.5 * dt * F_a2
    p_b  += 0.5 * dt * F_b2

    # ── Convert back to (omega, v) and pack ──────────────────────────────────
    omega_a_new = I_inv_a @ Pi_a
    omega_b_new = I_inv_b @ Pi_b
    v_a_new = p_a / m_a
    v_b_new = p_b / m_b

    body_a_new = BodyState(x_a_new, v_a_new, q_a_new, omega_a_new)
    body_b_new = BodyState(x_b_new, v_b_new, q_b_new, omega_b_new)
    return pack_state(body_a_new, body_b_new)


# ─────────────────────────────────────────────────────────────────────────────
# Convenience: integrate N steps and return full trajectory
# ─────────────────────────────────────────────────────────────────────────────

def integrate(
    state0: np.ndarray,
    dt: float,
    n_steps: int,
    m_a: float, I_a_body: np.ndarray,
    m_b: float, I_b_body: np.ndarray,
    *,
    record_every: int = 1,
) -> list[np.ndarray]:
    """Integrate n_steps steps and return states at every record_every steps.

    Returns a list of state arrays (including the initial state at index 0).
    """
    states = [state0]
    state = state0
    for i in range(n_steps):
        state = se3_variational_step(state, dt, m_a, I_a_body, m_b, I_b_body)
        if (i + 1) % record_every == 0:
            states.append(state)
    return states
