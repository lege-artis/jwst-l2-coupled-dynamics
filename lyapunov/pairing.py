"""
pairing.py — Hamiltonian-pairing gate helpers for the Phase 3 full spectrum.

Reuses lyapunov.reference.hamiltonian_pairing_check (DO NOT EDIT reference.py)
and adds two convenience gates:

  spectrum_sum_gate(exponents, threshold) — checks |Σλ| < threshold
  spurious_mode_absent(exponents, cutoff) — checks no exponent is < −cutoff

For the 12-dim deflated/Lie-alg spectrum of scenario (c):
  • sum_exponents must be ≪ 0.14 s⁻¹ (the R^14 contamination level)
  • The two spurious modes (≈ −0.064, −0.071 s⁻¹) must be absent
  • Hamiltonian pairing with rtol=0.5 must PASS for the inner pairs
  • n_near_zero ≥ 2 (energy + reduced-energy invariant directions)

PHASE 3 DELIVERABLE.

References:
  _audit/phase3-prep/FULL-SPECTRUM-FINDING-2026-06-07.md
  lyapunov/reference.py  (hamiltonian_pairing_check — DO NOT EDIT)
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

# Re-export the existing gate (no changes to reference.py).
from lyapunov.reference import hamiltonian_pairing_check  # noqa: F401


# ─────────────────────────────────────────────────────────────────────────────
# Phase 3 additional gates
# ─────────────────────────────────────────────────────────────────────────────

def spectrum_sum_gate(exponents: np.ndarray, threshold: float = 1e-2) -> bool:
    """AG-P3-1 primary: |Σλ| must be below threshold.

    For a Hamiltonian flow the true spectrum sums to zero.  The R^14 flat
    embedding gave Σλ = −0.14 s⁻¹ due to two spurious contracting modes.
    After deflation or in the Lie-alg embedding the sum should be ≪ 0.14,
    and for a well-converged run should satisfy |Σλ| < 1e-2.

    Args:
        exponents:  array of Lyapunov exponents (12-dim spectrum expected).
        threshold:  acceptance threshold (default 1e-2 s⁻¹).

    Returns True iff |sum(exponents)| < threshold.
    """
    return float(abs(np.sum(exponents))) < threshold


def spurious_modes_absent(
    exponents: np.ndarray,
    cutoff: float = 1e-2,
) -> bool:
    """AG-P3-1 secondary: no exponent more negative than −cutoff.

    The two spurious R^14 norm modes are at ≈ −0.064 and −0.071 s⁻¹.
    A clean 12-dim spectrum should have no exponent below −cutoff = −1e-2.
    The most-negative dynamical exponent in scenario (c) is expected near −λ_max
    (the Hamiltonian-pairing partner), which is of order 1e-5 to 1e-3.

    Args:
        exponents: array of Lyapunov exponents in descending order.
        cutoff:    magnitude threshold (default 1e-2 s⁻¹).

    Returns True iff all exponents > −cutoff.
    """
    return bool(np.all(exponents > -cutoff))


def full_spectrum_gates(
    exponents: np.ndarray,
    *,
    sum_threshold: float = 1e-2,
    spurious_cutoff: float = 1e-2,
    pairing_rtol: float = 0.5,
) -> dict:
    """Run all Phase 3 full-spectrum gates and return a summary dict.

    Args:
        exponents:       12-element Lyapunov spectrum in descending order.
        sum_threshold:   |Σλ| < sum_threshold   (AG-P3-1 primary)
        spurious_cutoff: min(λ) > −spurious_cutoff  (AG-P3-1 secondary)
        pairing_rtol:    relative tolerance for hamiltonian_pairing_check (AG-P3-2)

    Returns dict:
        sum_gate_pass:       bool (AG-P3-1 primary)
        spurious_gate_pass:  bool (AG-P3-1 secondary)
        pairing:             full output of hamiltonian_pairing_check (AG-P3-2)
        sum_exponents:       float
        min_exponent:        float
    """
    pairing = hamiltonian_pairing_check(exponents, rtol=pairing_rtol)
    return {
        'sum_gate_pass':      spectrum_sum_gate(exponents, sum_threshold),
        'spurious_gate_pass': spurious_modes_absent(exponents, spurious_cutoff),
        'pairing':            pairing,
        'sum_exponents':      float(np.sum(exponents)),
        'min_exponent':       float(np.min(exponents)),
    }
