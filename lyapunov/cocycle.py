"""
cocycle.py — Euclidean R^14 cocycle propagator for Lyapunov-spectrum computation.

OQ-PHASE2-2 DEFAULT = Euclidean: linearise the equations of motion in the
quaternion-component + translation-component R^14 chart and propagate tangent
vectors in R^14.  Reorthogonalise via Gram-Schmidt in R^14.

State-space: R^14 = q_a(4) + ω_a(3) + q_b(4) + ω_b(3)
  (translations omitted for purely rotational Lyapunov study; include all 26
  components for the full coupled translational+rotational spectrum)

SITTING 3 DELIVERABLE.

References:
  Benettin et al. (1980) Meccanica 15, 9–30.
  Skokos (2010) Lect. Notes Phys. 790.
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))


# ─────────────────────────────────────────────────────────────────────────────
# R^14 ↔ R^26 embedding / extraction helpers
# ─────────────────────────────────────────────────────────────────────────────
# R^26 state layout (one rigid body = 13 components):
#   Body A: x_a[0:3], v_a[3:6], q_a[6:10], ω_a[10:13]
#   Body B: x_b[13:16], v_b[16:19], q_b[19:23], ω_b[23:26]
#
# R^14 = q_a(4) + ω_a(3) + q_b(4) + ω_b(3)
#   R^14[0:4]   ↔ R^26[6:10]   (q_a)
#   R^14[4:7]   ↔ R^26[10:13]  (ω_a)
#   R^14[7:11]  ↔ R^26[19:23]  (q_b)
#   R^14[11:14] ↔ R^26[23:26]  (ω_b)

# (t_start, t_end, s_start, s_end) — tangent slice → state slice
_R14_MAP = [
    (0,  4,  6, 10),
    (4,  7, 10, 13),
    (7, 11, 19, 23),
    (11, 14, 23, 26),
]


def _embed_r14(v14: np.ndarray) -> np.ndarray:
    """Embed R^14 rotational tangent into R^26 (zero-fill x/v slots)."""
    v26 = np.zeros(26)
    for t0, t1, s0, s1 in _R14_MAP:
        v26[s0:s1] = v14[t0:t1]
    return v26


def _extract_r14(v26: np.ndarray) -> np.ndarray:
    """Extract R^14 rotational components from an R^26 vector."""
    return np.concatenate([v26[s0:s1] for _, _, s0, s1 in _R14_MAP])


# ─────────────────────────────────────────────────────────────────────────────
# Cocycle propagator
# ─────────────────────────────────────────────────────────────────────────────

def propagate_tangent(
    state: np.ndarray,
    tangent: np.ndarray,
    dt: float,
    m_a: float, I_a_body: np.ndarray,
    m_b: float, I_b_body: np.ndarray,
    *,
    fd_eps: float = 1e-7,
) -> np.ndarray:
    """Propagate one tangent vector by one integrator step via finite differences.

    Computes the Jacobian-vector product using the symmetric finite-difference
    approximation:
      Dφ(x) · v ≈ [φ(x + ε·v̂) − φ(x − ε·v̂)] / (2ε) · ‖v‖

    where v̂ = v / ‖v‖ is the unit tangent and φ is one KDK integrator step.
    Normalising before perturbation keeps the FD perturbation magnitude
    independent of the tangent-vector norm (which can grow during the BGGS run).

    tangent may be in R^14 (rotational subspace) or R^26 (full state).  The
    embedding / extraction is handled transparently.

    Note: the integrator normalises quaternions internally, so the slight
    quaternion off-normality introduced by the FD perturbation is corrected
    each step.

    Args:
        state:    26-component reference state at the current time t.
        tangent:  tangent vector (R^14 or R^26).
        dt:       integrator timestep (seconds).
        fd_eps:   finite-difference step size (default 1e-7; near-optimal for
                  double-precision states of O(1) magnitude).

    Returns:
        Propagated tangent in the same space as the input (R^14 or R^26).
    """
    from integrators.symplectic_se3_variational import se3_variational_step

    n = len(tangent)
    if n == 14:
        v_full = _embed_r14(tangent)
    elif n == 26:
        v_full = tangent.copy()
    else:
        raise ValueError(f"tangent must be R^14 or R^26, got dim={n}")

    # Normalise perturbation direction to avoid amplitude-dependent FD error.
    norm_v = np.linalg.norm(v_full)
    if norm_v < 1e-300:
        # Degenerate tangent — return zero vector; GS handles re-injection.
        return np.zeros(n)
    v_unit = v_full / norm_v

    # Symmetric FD: J · v_unit ≈ [φ(x + ε·v_unit) − φ(x − ε·v_unit)] / (2ε)
    phi_p = se3_variational_step(state + fd_eps * v_unit, dt,
                                 m_a, I_a_body, m_b, I_b_body)
    phi_m = se3_variational_step(state - fd_eps * v_unit, dt,
                                 m_a, I_a_body, m_b, I_b_body)
    d_full = (phi_p - phi_m) / (2.0 * fd_eps)

    # Restore original tangent magnitude: J · v = (J · v_unit) · ‖v‖
    d_full *= norm_v

    return _extract_r14(d_full) if n == 14 else d_full


# ─────────────────────────────────────────────────────────────────────────────
# Modified Gram-Schmidt orthonormalisation
# ─────────────────────────────────────────────────────────────────────────────

def gram_schmidt(basis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Modified Gram-Schmidt orthonormalisation of columns of basis ∈ R^{n×k}.

    Uses the modified (column-by-column) variant which is numerically more
    stable than the classical (row-by-row) variant.

    Returns:
        Q:          R^{n×k} array with orthonormal columns.
        log_norms:  R^k array where log_norms[i] = log(‖v_i‖ before
                    normalisation), i.e. the logarithm of the i-th stretching
                    factor.  Used by the BGGS algorithm to accumulate Lyapunov
                    exponents.  Returns −∞ for degenerate (zero-norm) columns,
                    which are replaced by a canonical basis vector.

    References:
        Benettin et al. (1980) Meccanica 15, §2 (GS reorthonormalisation step).
        Golub & Van Loan (2013) "Matrix Computations" §5.2 (modified GS).
    """
    n, k = basis.shape
    Q = basis.copy().astype(float)
    log_norms = np.zeros(k)

    for i in range(k):
        # Modified GS: subtract projections onto already-orthonormal columns.
        for j in range(i):
            Q[:, i] -= np.dot(Q[:, j], Q[:, i]) * Q[:, j]

        norm_i = np.linalg.norm(Q[:, i])

        if norm_i > 1e-300:
            log_norms[i] = np.log(norm_i)
            Q[:, i] /= norm_i
        else:
            # Degenerate column: replace with a canonical basis vector not
            # yet represented in the orthonormal set.
            log_norms[i] = -np.inf
            replacement = np.zeros(n)
            for candidate in range(n):
                e = np.zeros(n)
                e[candidate] = 1.0
                # Subtract projections onto existing orthonormal columns
                for j in range(i):
                    e -= np.dot(Q[:, j], e) * Q[:, j]
                norm_e = np.linalg.norm(e)
                if norm_e > 0.5:
                    replacement = e / norm_e
                    break
            Q[:, i] = replacement

    return Q, log_norms
