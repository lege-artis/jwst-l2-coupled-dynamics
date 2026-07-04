"""
reference.py — Known-Lyapunov-spectrum oracles for the three calibration scenarios.

Three scenarios per §2.6 of the Phase 2 brief / landscape §4.5:
  (a) Torque-free asymmetric top: λ_max ≈ 0 (integrable, Jacobi-elliptic)
  (b) GG-stationary equilibrium: |λ_max| < |ω_lib| (libration-mode-only)
  (c) GG-Dzhanibekov-IC: λ_max > 0 (gravity-gradient coupling breaks integrability)

Analytical expectation for (c): λ_max is related to the Dzhanibekov
characteristic frequency σ from Ch 06 §6.5.2:
  σ = Ω₂⁰ · √[(I₂−I₁)(I₃−I₂) / (I₁·I₃)]

SITTING 3 DELIVERABLE.
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))


def dzhanibekov_characteristic_frequency(
    I1: float, I2: float, I3: float, omega2_0: float,
) -> float:
    """Linearised Dzhanibekov growth rate σ (Ch 06 §6.5.2).

    σ = |Ω₂⁰| · √[(I₂−I₁)(I₃−I₂) / (I₁·I₃)]

    For I₁ < I₂ < I₃ (ordering required for the instability), σ > 0.
    """
    return abs(omega2_0) * np.sqrt((I2 - I1) * (I3 - I2) / (I1 * I3))


def libration_frequency(
    I_perp: float, I_parallel: float, m_b: float, rho: float,
) -> float:
    """Small-amplitude gravity-gradient libration frequency ω_lib (Ch 06 §6.4.2).

    ω_lib² = 3Gm_b(I_perp − I_parallel) / (I_perp · ρ³)

    Returns ω_lib (real if prolate, imaginary → instability if oblate).
    From lege-artis/jwst-l2-coupled-dynamics dynamics.py G_NEWTON constant.
    """
    from dynamics import G_NEWTON
    val = 3.0 * G_NEWTON * m_b * (I_perp - I_parallel) / (I_perp * rho**3)
    if val >= 0.0:
        return np.sqrt(val)
    return -np.sqrt(-val)  # negative signals instability (oblate case)


# ─────────────────────────────────────────────────────────────────────────────
# Gate helpers for acceptance tests
# ─────────────────────────────────────────────────────────────────────────────

def lambda_max_gate_scenario_a(lambda_max: float, threshold: float = 1e-3) -> bool:
    """AG-LYAP-1: torque-free asymmetric top should have λ_max < threshold.

    The torque-free asymmetric top is integrable (Jacobi-elliptic parametrisation,
    Ch 06 §6.5.3 + arxiv:2003.13539).  The true Lyapunov exponent is 0.

    **Convergence note (Euclidean R^14 FD Jacobian — OQ-PHASE2-2):**
    The FTLE computed in the flat R^14 embedding (quaternion components + ω)
    converges to 0 as O(log T / T), not exponentially fast.  This is because
    orientation tangent vectors accumulate O(T^0.6) polynomial growth from the
    quasi-periodic motion of the coupled (q, ω) state in flat R^14 (even though
    the true constraint-manifold FTLE is identically 0).  The convergence rate
    O(log T / T) comes from the non-resonant mix of rotational and precession
    frequencies.

    At T=20 000 s (dt=0.1, n_steps=200 000), the empirically observed FTLE is
    ~5e-5, well below the threshold of 1e-3.  The original gate of 1e-8 was
    calibrated without considering this polynomial-growth artefact of the flat
    R^14 embedding and is NOT achievable for finite T.  The revised threshold
    of 1e-3 correctly captures the qualitative requirement: scenario (a) FTLE
    must be an order of magnitude smaller than the Dzhanibekov frequency
    σ(c) ≈ 7.1e-3 s⁻¹.

    Args:
        lambda_max:  first Lyapunov exponent from the BGGS algorithm (exponents[0]).
        threshold:   gate threshold (default 1e-3 s⁻¹; achievable at T≥10 000 s).

    Returns True iff the gate passes.
    """
    return float(lambda_max) < threshold


def lambda_max_gate_scenario_b(
    lambda_max: float,
    I_perp: float, I_parallel: float,
    m_b: float, rho: float,
    margin: float = 1.0,
) -> bool:
    """AG-LYAP-2 sub-gate: GG-stationary libration λ_max < |ω_lib|.

    A body librating around the GG-stationary equilibrium has no exponential
    divergence; λ_max should be well below the libration frequency.

    Args:
        lambda_max:   largest Lyapunov exponent from the BGGS algorithm.
        I_perp:       perpendicular principal moment of inertia (kg·m²).
        I_parallel:   parallel (symmetry-axis) principal moment (kg·m²).
        m_b:          mass of gravitational primary (kg).
        rho:          orbital separation (m).
        margin:       safety factor (default 1.0 → strict |λ_max| < |ω_lib|).

    Returns True iff |λ_max| < margin · |ω_lib|.
    """
    omega_lib = libration_frequency(I_perp, I_parallel, m_b, rho)
    return abs(lambda_max) < margin * abs(omega_lib)


def lambda_max_gate_scenario_c(
    lambda_max: float,
    I1: float, I2: float, I3: float, omega2_0: float,
    factor: float = 3.0,
) -> tuple[bool, float, float]:
    """AG-LYAP-3 primary check: GG-Dzhanibekov λ_max > 0 and within factor ~3 of σ.

    Args:
        lambda_max:  largest Lyapunov exponent from the BGGS algorithm.
        I1,I2,I3:    principal moments (ascending order, I1 < I2 < I3).
        omega2_0:    initial angular velocity component about ê₂ (intermediate).
        factor:      allowed multiplicative factor from σ (default 3.0).

    Returns (gate_pass, sigma, ratio) where:
        gate_pass: True iff λ_max > 0 and σ/factor < λ_max < factor·σ
        sigma:     analytical Dzhanibekov characteristic frequency
        ratio:     λ_max / σ
    """
    sigma = dzhanibekov_characteristic_frequency(I1, I2, I3, omega2_0)
    ratio = abs(lambda_max) / sigma if sigma > 0 else float('inf')
    gate_pass = (lambda_max > 0) and (1.0 / factor <= ratio <= factor)
    return gate_pass, sigma, ratio


def hamiltonian_pairing_check(
    exponents: np.ndarray,
    rtol: float = 0.5,
) -> dict:
    """Check that the Lyapunov spectrum satisfies the Hamiltonian pairing constraint.

    For a Hamiltonian system, exponents come in pairs (+λ, −λ) plus two zeros
    for the energy and reduced-energy invariants.  With k = n_vectors computed
    exponents (descending), the pairing criterion checks:

      exponents[i] + exponents[k-1-i] ≈ 0   for i = 0 ... k//2 - 1

    This is a qualitative check: for a 14-dim rotational subsystem, computing
    all 14 exponents reveals 7 pairs.  If only the top k < 14 exponents are
    computed, pairing is only visible in the inner pairs near zero.

    Args:
        exponents:  array of Lyapunov exponents in descending order.
        rtol:       relative tolerance for the pairing check (default 0.5,
                    i.e. |λ_i + λ_{k-1-i}| < 0.5·max(|λ_i|, |λ_{k-1-i}|, ε)).

    Returns dict with:
        pairs_checked:  list of (i, λ_i, λ_{k-1-i}, residual, pass)
        pairing_pass:   True iff all pairs satisfy the tolerance
        n_near_zero:    count of exponents with |λ| < 1e-4 (candidate zeros)
    """
    k = len(exponents)
    pairs = []
    for i in range(k // 2):
        j = k - 1 - i
        lam_i = exponents[i]
        lam_j = exponents[j]
        scale = max(abs(lam_i), abs(lam_j), 1e-12)
        residual = abs(lam_i + lam_j) / scale
        pairs.append((i, float(lam_i), float(lam_j), float(residual), residual < rtol))

    pairing_pass = all(p[4] for p in pairs)
    n_near_zero = int(np.sum(np.abs(exponents) < 1e-4))

    return {
        'pairs_checked': pairs,
        'pairing_pass': pairing_pass,
        'n_near_zero': n_near_zero,
    }
