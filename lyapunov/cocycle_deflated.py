"""
cocycle_deflated.py — R^14 FD cocycle with quaternion-radial deflation (Tier A).

Phase 3 interim approach: keep the flat R^14 FD cocycle from cocycle.py but
deflate tangent vectors against the two quaternion-norm (radial) directions
  r̂_A = q_A / |q_A|  (in R^14 indices 0:4)
  r̂_B = q_B / |q_B|  (in R^14 indices 7:11)
before each Gram-Schmidt step, recovering a 12-dimensional effective spectrum
without a full Lie-algebraic rewrite.

This eliminates the two spurious large-negative modes (≈ -0.064, -0.071 s⁻¹)
identified in _audit/phase3-prep/FULL-SPECTRUM-FINDING-2026-06-07.md, restoring
the Hamiltonian pairing structure visible in the 12-dim subspace.

State-space (R^14 layout, inherited from cocycle.py):
  R^14[0:4]   = q_A   ← deflation target for body A
  R^14[4:7]   = ω_A
  R^14[7:11]  = q_B   ← deflation target for body B
  R^14[11:14] = ω_B

PHASE 3 TIER A DELIVERABLE.

References:
  _audit/phase3-prep/FULL-SPECTRUM-FINDING-2026-06-07.md  (empirical motivation)
  lyapunov/cocycle.py  (R^14 FD cocycle baseline — DO NOT EDIT)
  Benettin et al. (1980) Meccanica 15, 9-30
  Skokos, C. (2010) Lect. Notes Phys. 790, §3.1
"""
from __future__ import annotations
import sys
from pathlib import Path
from typing import Callable

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from lyapunov.cocycle import propagate_tangent, gram_schmidt


# ─────────────────────────────────────────────────────────────────────────────
# Quaternion-radial deflation helpers
# ─────────────────────────────────────────────────────────────────────────────
# R^26 state layout (from cocycle.py):
#   q_A at [6:10],  ω_A at [10:13]
#   q_B at [19:23], ω_B at [23:26]
# R^14 quaternion slots:
#   R^14[0:4]  ↔ q_A
#   R^14[7:11] ↔ q_B


def _radial_dirs(state: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Unit quaternion-radial directions for body A and B, embedded in R^14.

    Returns (r_A, r_B), each a unit vector in R^14 non-zero only in the
    quaternion slot of the respective body.  The CURRENT state's quaternions
    are used so that the radial directions co-evolve with the trajectory.
    """
    q_a = state[6:10]
    q_b = state[19:23]

    r_A = np.zeros(14)
    r_A[0:4] = q_a / np.linalg.norm(q_a)

    r_B = np.zeros(14)
    r_B[7:11] = q_b / np.linalg.norm(q_b)

    return r_A, r_B


def deflate_tangent(
    v: np.ndarray,
    r_A: np.ndarray,
    r_B: np.ndarray,
) -> np.ndarray:
    """Project R^14 tangent vector v orthogonal to the two radial directions.

    Removes the quaternion-norm components that point off the S^3 constraint
    manifold, leaving v in the 12-dimensional dynamical subspace.
    """
    v = v - np.dot(v, r_A) * r_A
    v = v - np.dot(v, r_B) * r_B
    return v


# ─────────────────────────────────────────────────────────────────────────────
# BGGS driver — deflated cocycle
# ─────────────────────────────────────────────────────────────────────────────

def lyapunov_spectrum_deflated(
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
) -> dict:
    """BGGS Lyapunov spectrum using the R^14 FD cocycle + quaternion-radial deflation.

    Runs the standard BGGS algorithm in R^14 but deflates each tangent vector
    against the two quaternion-norm directions (r̂_A, r̂_B) after every propagation
    step and before Gram-Schmidt, restricting the computation to the 12-dimensional
    dynamical subspace.  The deflation directions co-evolve with the trajectory.

    Args:
        state0:    initial 26-component state
        step_fn:   integrator: step_fn(state, dt, m_a, I_a_body, m_b, I_b_body)
        dt:        timestep (seconds)
        n_steps:   total integration steps for the production run
        n_vectors: number of tangent vectors; must be ≤ 12
        tau_qr:    QR-cadence (None → autotuned from pilot run)

    Returns dict with keys:
        exponents:         n_vectors Lyapunov exponents in descending order
        sum_exponents:     sum of all computed exponents (should be ≪ 0.14 for 12-dim)
        tau_qr_used:       production QR cadence (s)
        lambda_max_pilot:  λ_max estimate from pilot run
        n_qr:              number of QR reorthogonalisations
        running_exponents: list of (t, exponents) snapshots at each QR step
        T_total:           total integration time in production run (s)
        final_state:       26-component state after production run
    """
    if I_a_body is None:
        I_a_body = np.eye(3)
    if I_b_body is None:
        I_b_body = np.eye(3)
    if n_vectors > 12:
        raise ValueError(
            f"n_vectors={n_vectors} exceeds the 12-dim deflated subspace"
        )

    # ── Initialise basis in the deflated 12-dim subspace ─────────────────────
    rng = np.random.default_rng(seed=42)
    raw = rng.standard_normal((14, n_vectors))
    r_A0, r_B0 = _radial_dirs(state0)
    for j in range(n_vectors):
        raw[:, j] = deflate_tangent(raw[:, j], r_A0, r_B0)
    Q, _ = gram_schmidt(raw)

    state = state0.copy()

    # ── Pilot run: estimate λ_max ─────────────────────────────────────────────
    if tau_qr is None:
        tau_qr_pilot = 10.0 * dt
        steps_pilot = max(1, int(round(tau_qr_pilot / dt)))
        n_pilot_epochs = 20

        accum_pilot = np.zeros(n_vectors)
        T_pilot = 0.0

        for _ in range(n_pilot_epochs):
            for _ in range(steps_pilot):
                r_A, r_B = _radial_dirs(state)
                for j in range(n_vectors):
                    Qj = propagate_tangent(
                        state, Q[:, j], dt, m_a, I_a_body, m_b, I_b_body,
                    )
                    Q[:, j] = deflate_tangent(Qj, r_A, r_B)
                state = step_fn(state, dt, m_a, I_a_body, m_b, I_b_body)
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
            r_A, r_B = _radial_dirs(state)
            for j in range(n_vectors):
                Qj = propagate_tangent(
                    state, Q[:, j], dt, m_a, I_a_body, m_b, I_b_body,
                )
                Q[:, j] = deflate_tangent(Qj, r_A, r_B)
            state = step_fn(state, dt, m_a, I_a_body, m_b, I_b_body)
            T_total += dt
        Q, log_norms = gram_schmidt(Q)
        accum += log_norms
        n_qr += 1
        if T_total > 0:
            running_exponents.append((T_total, accum / T_total))

    if remainder > 0:
        for _ in range(remainder):
            r_A, r_B = _radial_dirs(state)
            for j in range(n_vectors):
                Qj = propagate_tangent(
                    state, Q[:, j], dt, m_a, I_a_body, m_b, I_b_body,
                )
                Q[:, j] = deflate_tangent(Qj, r_A, r_B)
            state = step_fn(state, dt, m_a, I_a_body, m_b, I_b_body)
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
