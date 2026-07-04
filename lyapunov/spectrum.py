"""
spectrum.py — Lyapunov spectrum via Benettin-Galgani-Giorgilli-Strelcyn 1980.

Algorithm:
  1. Run a pilot integration (tau_QR_pilot = 10 · dt) to estimate λ_max.
  2. Set production QR cadence: tau_QR = 0.1 / λ_max_pilot (clamped).
  3. Production run: propagate state + k tangent vectors; reorthogonalise
     (Gram-Schmidt) every tau_QR; accumulate log-norms.
  4. Divide accumulated log-norms by total integration time → spectrum.

Phase 2 gates (landscape §4.4):
  • Qualitative: sign of λ_max per scenario
  • Ordering: Hamiltonian pairing (each +λ paired with −λ, two zeros)
  • Order-of-magnitude: λ_max within factor ~3 of analytical expectation

SITTING 3 DELIVERABLE.

References:
  Benettin, G., Galgani, L., Giorgilli, A., Strelcyn, J.-M. (1980)
    Meccanica 15, 9–30.
  Skokos, C. (2010) Lect. Notes Phys. 790, Springer.
  Geist, K., Parlitz, U., Lauterborn, W. (1990)
    Prog. Theor. Phys. 83, 875–893.
"""
from __future__ import annotations
import sys
from pathlib import Path
from typing import Callable

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))


# ─────────────────────────────────────────────────────────────────────────────
# Main BGGS driver
# ─────────────────────────────────────────────────────────────────────────────

def lyapunov_spectrum(
    state0: np.ndarray,
    step_fn: Callable,
    dt: float,
    n_steps: int,
    n_vectors: int,
    m_a: float, I_a_body: np.ndarray,
    m_b: float, I_b_body: np.ndarray,
    *,
    tau_qr_pilot: float | None = None,
    tau_qr: float | None = None,
) -> dict:
    """Compute the Lyapunov spectrum for the two-body coupled dynamics.

    Implements the Benettin-Galgani-Giorgilli-Strelcyn (1980) algorithm with
    pilot-then-production QR-cadence autotuning (landscape §4.3).

    Tangent-space selection:
      n_vectors ≤ 14 → work in R^14 (rotational subspace; OQ-PHASE2-2 default)
      n_vectors > 14 → work in R^26 (full state)

    The pilot run (n_pilot_epochs × tau_qr_pilot per epoch) estimates λ_max.
    The production QR cadence is set to 0.1 / |λ_max_pilot| (clamped to
    [10·dt, n_steps·dt/10] to handle near-zero pilot estimates).

    Args:
        state0:        initial 26-component state
        step_fn:       integrator step function, called as
                         step_fn(state, dt, m_a, I_a_body, m_b, I_b_body)
        dt:            timestep (seconds)
        n_steps:       total integration steps for production run
        n_vectors:     number of tangent vectors (≤ 26; use ≤ 14 for rotational)
        tau_qr_pilot:  QR-cadence for pilot run (None → 10·dt)
        tau_qr:        QR-cadence for production run (None → autotuned)

    Returns dict with keys:
        exponents:         R^{n_vectors} Lyapunov exponents in descending order
        tau_qr_used:       production QR cadence (seconds)
        lambda_max_pilot:  λ_max estimate from pilot run
        n_qr:              number of QR reorthogonalisations performed
        running_exponents: list of (t, exponents) snapshots at each QR step
        T_total:           total integration time in production run (seconds)
        final_state:       26-component state after production run
    """
    from lyapunov.cocycle import propagate_tangent, gram_schmidt

    # ── Dimension selection ──────────────────────────────────────────────────
    dim = 14 if n_vectors <= 14 else 26
    if n_vectors > dim:
        raise ValueError(f"n_vectors={n_vectors} exceeds tangent dim={dim}")

    # ── Initialise random orthonormal tangent basis ──────────────────────────
    rng = np.random.default_rng(seed=42)
    raw = rng.standard_normal((dim, n_vectors))
    Q, _ = gram_schmidt(raw)

    state = state0.copy()

    # ─── Pilot run: estimate λ_max ──────────────────────────────────────────
    if tau_qr is None:
        if tau_qr_pilot is None:
            tau_qr_pilot = 10.0 * dt

        steps_per_qr_pilot = max(1, int(round(tau_qr_pilot / dt)))
        n_pilot_epochs = 20  # lightweight: 20 × 10 = 200 steps total

        accum_pilot = np.zeros(n_vectors)
        T_pilot = 0.0

        for _ in range(n_pilot_epochs):
            # Propagate state + tangent vectors for one QR epoch
            for _ in range(steps_per_qr_pilot):
                for j in range(n_vectors):
                    Q[:, j] = propagate_tangent(
                        state, Q[:, j], dt,
                        m_a, I_a_body, m_b, I_b_body,
                    )
                state = step_fn(state, dt, m_a, I_a_body, m_b, I_b_body)
                T_pilot += dt
            # Gram-Schmidt + accumulate
            Q, log_norms = gram_schmidt(Q)
            accum_pilot += log_norms

        # λ_max_pilot = total accumulated log-norm / total time
        lambda_max_pilot = float(accum_pilot[0] / T_pilot) if T_pilot > 0 else 0.0

        # Production QR cadence: target ~10 oscillations per epoch.
        # Clamp to [10·dt, n_steps·dt / 10] to avoid degenerate cases.
        t_total_prod = n_steps * dt
        if abs(lambda_max_pilot) > 1e-12:
            tau_qr = 0.1 / abs(lambda_max_pilot)
        else:
            # λ_max ≈ 0 (integrable): use a generous cadence
            tau_qr = min(100.0 * dt, t_total_prod / 5.0)

        tau_qr = max(10.0 * dt, min(tau_qr, t_total_prod / 3.0))
    else:
        lambda_max_pilot = 0.0

    # ─── Production run ──────────────────────────────────────────────────────
    steps_per_qr = max(1, int(round(tau_qr / dt)))
    accum = np.zeros(n_vectors)
    T_total = 0.0
    n_qr = 0
    running_exponents: list[tuple[float, np.ndarray]] = []

    n_complete = n_steps // steps_per_qr
    remainder = n_steps % steps_per_qr

    for _ in range(n_complete):
        for _ in range(steps_per_qr):
            for j in range(n_vectors):
                Q[:, j] = propagate_tangent(
                    state, Q[:, j], dt,
                    m_a, I_a_body, m_b, I_b_body,
                )
            state = step_fn(state, dt, m_a, I_a_body, m_b, I_b_body)
            T_total += dt
        Q, log_norms = gram_schmidt(Q)
        accum += log_norms
        n_qr += 1
        if T_total > 0:
            running_exponents.append((T_total, accum / T_total))

    # Handle remaining steps (partial epoch)
    if remainder > 0:
        for _ in range(remainder):
            for j in range(n_vectors):
                Q[:, j] = propagate_tangent(
                    state, Q[:, j], dt,
                    m_a, I_a_body, m_b, I_b_body,
                )
            state = step_fn(state, dt, m_a, I_a_body, m_b, I_b_body)
            T_total += dt
        # Final partial-epoch QR
        Q, log_norms = gram_schmidt(Q)
        accum += log_norms
        n_qr += 1
        if T_total > 0:
            running_exponents.append((T_total, accum / T_total))

    exponents = accum / T_total if T_total > 0 else accum.copy()

    return {
        'exponents': exponents,
        'tau_qr_used': tau_qr,
        'lambda_max_pilot': lambda_max_pilot,
        'n_qr': n_qr,
        'running_exponents': running_exponents,
        'T_total': T_total,
        'final_state': state,
    }
