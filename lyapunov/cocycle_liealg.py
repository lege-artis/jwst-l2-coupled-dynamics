"""
cocycle_liealg.py — Lie-algebraic so(3)×R^3 cocycle (Tier B).

Phase 3 canonical fix: replace the flat R^14 FD Jacobian with a tangent
propagator that works in the body-frame Lie algebra
  so(3) × R^3 per body  →  12-dimensional by construction.

Tangent state layout (R^12):
  [0:3]  ξ_A  ∈ so(3) ≅ R^3  — rotation-vector perturbation for body A
  [3:6]  δω_A ∈ R^3           — angular-velocity perturbation for body A
  [6:9]  ξ_B  ∈ so(3) ≅ R^3
  [9:12] δω_B ∈ R^3

Convention — RIGHT trivialisation (matches integrators/symplectic_se3_variational.py):
  orientation perturbation:  q_pert = q_normalise(q ⊗ rotvec_to_quat(ε · ξ̂))
  angular-velocity pert:     ω_pert = ω + ε · δω̂
  log map (result → so(3)): ξ_out  = 2·arctan2(|q_vec|, q_w) · q_vec/|q_vec|
                              where q_out = q_conj(q_ref) ⊗ q_step_result

The right-trivialisation is verified from symplectic_se3_variational.py line
  q_new = q_normalise(quat_multiply(q, rotvec_to_quat(v)))
where v = dt · ω_mid is the RIGHT-applied rotation vector.  The log map is
therefore  ξ_out = log(q_ref⁻¹ ⊗ q_pert_out) with the conjugate on the LEFT.

PHASE 3 TIER B DELIVERABLE.

References:
  Hairer, Lubich, Wanner (2006) Geometric Numerical Integration, §IV (Lie-group methods)
  Marsden, Ratiu (1999) Introduction to Mechanics and Symmetry, §9 (so(3) reduction)
  integrators/symplectic_se3_variational.py (convention anchor)
  lyapunov/cocycle.py (R^14 FD baseline — DO NOT EDIT)
  Benettin et al. (1980) Meccanica 15, 9-30
"""
from __future__ import annotations
import sys
from pathlib import Path
from typing import Callable

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from lyapunov.cocycle import gram_schmidt  # reuse; DO NOT EDIT cocycle.py


# ─────────────────────────────────────────────────────────────────────────────
# Quaternion helpers (right-trivialised convention)
# ─────────────────────────────────────────────────────────────────────────────

def _quat_conj(q: np.ndarray) -> np.ndarray:
    """Conjugate of unit quaternion (qw, qx, qy, qz)."""
    return np.array([q[0], -q[1], -q[2], -q[3]])


def _quat_mul(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Hamilton product p ⊗ q  (qw, qx, qy, qz convention)."""
    pw, px, py, pz = p
    qw, qx, qy, qz = q
    return np.array([
        pw*qw - px*qx - py*qy - pz*qz,
        pw*qx + px*qw + py*qz - pz*qy,
        pw*qy - px*qz + py*qw + pz*qx,
        pw*qz + px*qy - py*qx + pz*qw,
    ])


def _rotvec_to_quat(v: np.ndarray) -> np.ndarray:
    """Unit quaternion (qw, qx, qy, qz) encoding rotation vector v = θ n̂."""
    theta = np.linalg.norm(v)
    if theta < 1e-14:
        q = np.array([1.0, 0.5 * v[0], 0.5 * v[1], 0.5 * v[2]])
    else:
        n = v / theta
        h = 0.5 * theta
        s = np.sin(h)
        q = np.array([np.cos(h), s*n[0], s*n[1], s*n[2]])
    return q / np.linalg.norm(q)


def _quat_log(q: np.ndarray) -> np.ndarray:
    """Rotation vector (R^3) corresponding to unit quaternion (qw, qx, qy, qz).

    Returns the rotation vector v such that rotvec_to_quat(v) ≈ q.
    For |q_vec| < 1e-14 returns 2·q_vec (first-order; q_w ≈ 1 so arctan2 ≈ |q_vec|).
    """
    qw = float(q[0])
    qvec = q[1:]
    norm_v = np.linalg.norm(qvec)
    if norm_v < 1e-14:
        return 2.0 * qvec
    return 2.0 * np.arctan2(norm_v, qw) * (qvec / norm_v)


# ─────────────────────────────────────────────────────────────────────────────
# State perturbation helpers (embed R^12 tangent → R^26 perturbation)
# ─────────────────────────────────────────────────────────────────────────────

# R^26 state layout:
#   body A:  x_a[0:3], v_a[3:6], q_a[6:10], ω_a[10:13]
#   body B:  x_b[13:16], v_b[16:19], q_b[19:23], ω_b[23:26]

def _perturb_state(state: np.ndarray, tangent: np.ndarray, eps: float) -> np.ndarray:
    """Return state perturbed by eps · tangent in the so(3)×R^3 chart.

    tangent ∈ R^12: [ξ_A(3), δω_A(3), ξ_B(3), δω_B(3)].
    Orientation perturbation uses the RIGHT exponential map: q_pert = q ⊗ exp(ε·ξ).
    ω perturbation is additive.
    Translational components (x, v) are unchanged (purely rotational cocycle).
    """
    s = state.copy()

    xi_A  = tangent[0:3]
    dw_A  = tangent[3:6]
    xi_B  = tangent[6:9]
    dw_B  = tangent[9:12]

    # Body A
    q_a = state[6:10]
    s[6:10]  = _quat_mul(q_a, _rotvec_to_quat(eps * xi_A))
    s[6:10] /= np.linalg.norm(s[6:10])
    s[10:13] = state[10:13] + eps * dw_A

    # Body B
    q_b = state[19:23]
    s[19:23]  = _quat_mul(q_b, _rotvec_to_quat(eps * xi_B))
    s[19:23] /= np.linalg.norm(s[19:23])
    s[23:26] = state[23:26] + eps * dw_B

    return s


def _extract_liealg(
    phi_p: np.ndarray,
    phi_m: np.ndarray,
    q_a_ref: np.ndarray,
    q_b_ref: np.ndarray,
    eps: float,
    norm_v: float,
) -> np.ndarray:
    """Map the symmetric-FD result back to R^12 so(3)×R^3 tangent.

    Args:
        phi_p, phi_m:  R^26 states from +ε and −ε perturbed steps.
        q_a_ref, q_b_ref: reference quaternions after the unperturbed step.
        eps:  finite-difference step size used for the perturbation.
        norm_v: magnitude of the original tangent vector (restored at the end).

    Returns R^12 tangent vector.
    """
    scale = norm_v / (2.0 * eps)

    # Body A — quaternion part: log(q_ref⁻¹ ⊗ q_p) − log(q_ref⁻¹ ⊗ q_m)
    q_a_p = phi_p[6:10]
    q_a_m = phi_m[6:10]
    dq_a_p = _quat_mul(_quat_conj(q_a_ref), q_a_p)
    dq_a_m = _quat_mul(_quat_conj(q_a_ref), q_a_m)
    xi_A_out = (_quat_log(dq_a_p) - _quat_log(dq_a_m)) * scale

    # Body A — angular-velocity part: plain FD
    dw_A_out = (phi_p[10:13] - phi_m[10:13]) * scale

    # Body B — quaternion part
    q_b_p = phi_p[19:23]
    q_b_m = phi_m[19:23]
    dq_b_p = _quat_mul(_quat_conj(q_b_ref), q_b_p)
    dq_b_m = _quat_mul(_quat_conj(q_b_ref), q_b_m)
    xi_B_out = (_quat_log(dq_b_p) - _quat_log(dq_b_m)) * scale

    # Body B — angular-velocity part
    dw_B_out = (phi_p[23:26] - phi_m[23:26]) * scale

    return np.concatenate([xi_A_out, dw_A_out, xi_B_out, dw_B_out])


# ─────────────────────────────────────────────────────────────────────────────
# Lie-algebraic tangent propagator
# ─────────────────────────────────────────────────────────────────────────────

def propagate_tangent_liealg(
    state: np.ndarray,
    state_next: np.ndarray,
    tangent: np.ndarray,
    dt: float,
    m_a: float,
    I_a_body: np.ndarray,
    m_b: float,
    I_b_body: np.ndarray,
    step_fn: Callable,
    *,
    fd_eps: float = 1e-7,
) -> np.ndarray:
    """Propagate one R^12 Lie-algebraic tangent vector by one integrator step.

    Uses symmetric finite differences in the so(3)×R^3 chart:
      tangent_out ≈ [D_Lie φ(x)] · tangent
    where the perturbation is applied via the RIGHT exponential map and the
    result is extracted via the left-trivialised log  ξ_out = log(q_ref⁻¹ ⊗ q_out).

    state_next must be the unperturbed reference step  φ(state, dt);  it is
    computed once per integration step by the BGGS driver and shared across all
    n_vectors tangent vectors, so the total evaluation count remains
    2·n_vectors + 1 per integration step.

    Args:
        state:      R^26 current state at time t.
        state_next: R^26 reference state at time t + dt (one unperturbed step).
        tangent:    R^12 tangent vector in the so(3)×R^3 chart.
        fd_eps:     finite-difference step size (default 1e-7).

    Returns R^12 propagated tangent.
    """
    norm_v = np.linalg.norm(tangent)
    if norm_v < 1e-300:
        return np.zeros(12)

    v_hat = tangent / norm_v

    state_p = _perturb_state(state, v_hat,  fd_eps)
    state_m = _perturb_state(state, v_hat, -fd_eps)

    phi_p = step_fn(state_p, dt, m_a, I_a_body, m_b, I_b_body)
    phi_m = step_fn(state_m, dt, m_a, I_a_body, m_b, I_b_body)

    q_a_ref = state_next[6:10]
    q_b_ref = state_next[19:23]

    return _extract_liealg(phi_p, phi_m, q_a_ref, q_b_ref, fd_eps, norm_v)


# ─────────────────────────────────────────────────────────────────────────────
# BGGS driver — Lie-algebraic cocycle
# ─────────────────────────────────────────────────────────────────────────────

def lyapunov_spectrum_liealg(
    state0: np.ndarray,
    step_fn: Callable,
    dt: float,
    n_steps: int,
    n_vectors: int = 12,
    m_a: float = 1.0,
    I_a_body: np.ndarray | None = None,
    m_b: float = 1.0,
    I_b_body: np.ndarray | None = None,
    *,
    tau_qr: float | None = None,
    fd_eps: float = 1e-7,
) -> dict:
    """BGGS Lyapunov spectrum using the so(3)×R^3 Lie-algebraic cocycle.

    The tangent space is 12-dimensional by construction — orientation perturbations
    live in so(3) and are applied via the right exponential map, so the quaternion-
    norm direction never appears.  The Hamiltonian-pairing structure (6 pairs of
    ±λ plus near-zero invariant directions) is therefore observable without
    post-hoc deflation.

    Each integration step calls step_fn once for the reference state and twice
    for each tangent vector (symmetric FD), for a total of 1 + 2·n_vectors
    evaluations per step — the same count as the R^14 flat cocycle.

    Args:
        state0:    initial 26-component state
        step_fn:   integrator: step_fn(state, dt, m_a, I_a_body, m_b, I_b_body)
        dt:        timestep (seconds)
        n_steps:   total integration steps
        n_vectors: number of tangent vectors; must be ≤ 12
        tau_qr:    QR-cadence (None → autotuned from pilot run)
        fd_eps:    finite-difference step size

    Returns dict with keys:
        exponents:         n_vectors Lyapunov exponents in descending order
        sum_exponents:     sum of exponents (should be ~0 for Hamiltonian flow)
        tau_qr_used, lambda_max_pilot, n_qr, running_exponents, T_total,
        final_state        (same as lyapunov_spectrum)
    """
    if I_a_body is None:
        I_a_body = np.eye(3)
    if I_b_body is None:
        I_b_body = np.eye(3)
    if n_vectors > 12:
        raise ValueError(
            f"n_vectors={n_vectors} exceeds the 12-dim Lie-algebraic tangent space"
        )

    # ── Initialise random orthonormal basis in R^12 ───────────────────────────
    rng = np.random.default_rng(seed=42)
    raw = rng.standard_normal((12, n_vectors))
    Q, _ = gram_schmidt(raw)

    state = state0.copy()

    _prop = lambda st, st_next, tv: propagate_tangent_liealg(
        st, st_next, tv, dt, m_a, I_a_body, m_b, I_b_body, step_fn,
        fd_eps=fd_eps,
    )

    # ── Pilot run ─────────────────────────────────────────────────────────────
    if tau_qr is None:
        tau_qr_pilot = 10.0 * dt
        steps_pilot = max(1, int(round(tau_qr_pilot / dt)))
        n_pilot_epochs = 20

        accum_pilot = np.zeros(n_vectors)
        T_pilot = 0.0

        for _ in range(n_pilot_epochs):
            for _ in range(steps_pilot):
                state_next = step_fn(state, dt, m_a, I_a_body, m_b, I_b_body)
                for j in range(n_vectors):
                    Q[:, j] = _prop(state, state_next, Q[:, j])
                state = state_next
                T_pilot += dt
            Q, log_norms = gram_schmidt(Q)
            accum_pilot += log_norms

        lambda_max_pilot = float(accum_pilot[0] / T_pilot) if T_pilot > 0 else 0.0
        t_prod = n_steps * dt
        if abs(lambda_max_pilot) > 1e-12:
            tau_qr = 0.1 / abs(lambda_max_pilot)
        else:
            tau_qr = min(100.0 * dt, t_prod / 5.0)
        tau_qr = max(10.0 * dt, min(tau_qr, t_prod / 3.0))
    else:
        lambda_max_pilot = 0.0

    # ── Production run ────────────────────────────────────────────────────────
    steps_per_qr = max(1, int(round(tau_qr / dt)))
    accum = np.zeros(n_vectors)
    T_total = 0.0
    n_qr = 0
    running_exponents: list[tuple[float, np.ndarray]] = []

    n_complete = n_steps // steps_per_qr
    remainder  = n_steps %  steps_per_qr

    for _ in range(n_complete):
        for _ in range(steps_per_qr):
            state_next = step_fn(state, dt, m_a, I_a_body, m_b, I_b_body)
            for j in range(n_vectors):
                Q[:, j] = _prop(state, state_next, Q[:, j])
            state = state_next
            T_total += dt
        Q, log_norms = gram_schmidt(Q)
        accum += log_norms
        n_qr += 1
        if T_total > 0:
            running_exponents.append((T_total, accum / T_total))

    if remainder > 0:
        for _ in range(remainder):
            state_next = step_fn(state, dt, m_a, I_a_body, m_b, I_b_body)
            for j in range(n_vectors):
                Q[:, j] = _prop(state, state_next, Q[:, j])
            state = state_next
            T_total += dt
        Q, log_norms = gram_schmidt(Q)
        accum += log_norms
        n_qr += 1
        if T_total > 0:
            running_exponents.append((T_total, accum / T_total))

    exponents = accum / T_total if T_total > 0 else accum.copy()

    return {
        'exponents': exponents,
        'sum_exponents': float(np.sum(exponents)),
        'tau_qr_used': tau_qr,
        'lambda_max_pilot': lambda_max_pilot,
        'n_qr': n_qr,
        'running_exponents': running_exponents,
        'T_total': T_total,
        'final_state': state,
    }
