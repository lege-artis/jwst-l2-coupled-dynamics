"""
symplectic_se3_splitting.py — Yoshida 4th-order splitting backstop.

OQ-PHASE2-1 DEFAULT = BOTH: this module is the Yoshida-splitting backstop
alongside symplectic_se3_variational.py.  SITTING 2 DELIVERABLE.

Algorithm: 4th-order Yoshida composition of two exactly-solvable flows
─────────────────────────────────────────────────────────────────────────────
The Hamiltonian is split H = T(p, Π) + V(x, q) where:

  T-flow (kinetic, drift substep for time s):
    • p and Π constant (momenta are integrals of the T-flow)
    • x → x + s·(p/m)                   (centre-of-mass drift)
    • Π → F^T·Π,  q → q ⊗ q(F)         (free rotation + coadjoint transport)
      where F = exp(s·ω),  ω = I⁻¹·Π   (body-frame free precession)
    The T-flow is EXACT: Π is constant so ω = I⁻¹·Π is constant, and
    exp(s·ω) integrates the free Euler equations exactly for time s.

  V-flow (gravitational, kick substep for time s):
    • x and q constant (the potential V does not depend on p or Π)
    • Π → Π + s·τ                        (gravity-gradient torque impulse)
    • p → p + s·F_grav                   (gravitational force impulse)

Yoshida 4th-order composition (Yoshida 1990 Table 1, n=2):
  w₁ = 1/(2−2^{1/3}),  w₀ = −2^{1/3}·w₁
  Composition: [drift c₁·dt] [kick d₁·dt] [drift c₂·dt] [kick d₂·dt]
               [drift c₃·dt] [kick d₃·dt] [drift c₄·dt]
  where:
    c₁ = c₄ = w₁/2
    c₂ = c₃ = (w₀+w₁)/2    (note: w₀ < 0 so c₂ < c₁)
    d₁ = d₃ = w₁
    d₂       = w₀           (negative weight — the Yoshida "backward" stage)

Properties:
  • 4th-order accurate (error ~ O(dt⁵) per HLW 2006 Ch VII)
  • Symplectic: each substep is symplectic; composition of symplectics is symplectic
  • Momentum-preserving: inertial L = R·Π preserved EXACTLY per drift substep
    (coadjoint transport preserves R·Π; V-kick doesn't change q/R)
  • Energy bounded oscillation (backward error analysis, HLW Ch IX)
  • SO(3) group structure preserved at every substep (no reprojection needed)
  • Same call signature as se3_variational_step — drop-in alternative

References:
  Yoshida, H. (1990) "Construction of high order symplectic integrators",
    Phys. Lett. A 150, 262–268.  [PRIMARY ANCHOR — coefficient derivation]
  Hairer, E., Lubich, C., Wanner, G. (2006) Geometric Numerical Integration,
    2nd ed., Springer, Ch VII.  [Splitting theory + Yoshida composition theory]
  Lee, T., Leok, M., McClamroch, N.H. (2007) CMDA 98, 121–144.  [Context]
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from dynamics import BodyState, pack_state, unpack_state, q_normalise

from integrators.symplectic_se3_variational import (
    rodrigues, rotvec_to_quat, quat_multiply, _forces_torques,
)


# ─────────────────────────────────────────────────────────────────────────────
# Yoshida 4th-order coefficients (Yoshida 1990 Table 1, n=2 case)
# ─────────────────────────────────────────────────────────────────────────────

_CBRT2 = 2.0 ** (1.0 / 3.0)
_W0 = -_CBRT2 / (2.0 - _CBRT2)        # ≈ −1.7024
_W1 = 1.0 / (2.0 - _CBRT2)             # ≈  1.3512

# Position-step (drift) coefficients: c[0..3]
YOSHIDA_C = np.array([
    _W1 / 2.0,
    (_W0 + _W1) / 2.0,
    (_W0 + _W1) / 2.0,
    _W1 / 2.0,
])

# Momentum-step (kick) coefficients: d[0..2]
YOSHIDA_D = np.array([_W1, _W0, _W1])


# ─────────────────────────────────────────────────────────────────────────────
# Substep helpers
# ─────────────────────────────────────────────────────────────────────────────

def _drift_substep(
    Pi: np.ndarray, q: np.ndarray, x: np.ndarray, p: np.ndarray,
    I_inv: np.ndarray, m: float, s: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Exact kinetic (T) flow for one body for time s.

    Π is constant ⟹ ω = I⁻¹·Π is constant ⟹ free rotation is exactly
    exp(s·ω).  Coadjoint transport Π_new = F^T·Π preserves R·Π exactly.

    Returns (Pi_new, q_new, x_new).
    """
    omega = I_inv @ Pi
    v = s * omega                              # rotation vector = angle·axis
    F = rodrigues(v)
    Pi_new = F.T @ Pi                          # coadjoint transport
    q_new = q_normalise(quat_multiply(q, rotvec_to_quat(v)))
    x_new = x + s * (p / m)
    return Pi_new, q_new, x_new


def _kick_substep(
    Pi_a: np.ndarray, Pi_b: np.ndarray,
    p_a: np.ndarray, p_b: np.ndarray,
    x_a: np.ndarray, x_b: np.ndarray,
    q_a: np.ndarray, q_b: np.ndarray,
    I_a_body: np.ndarray, I_b_body: np.ndarray,
    m_a: float, m_b: float,
    s: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Gravitational (V) flow substep for time s.

    x and q are frozen; only momenta change:
      Π → Π + s·τ    (torque impulse)
      p → p + s·F    (force impulse)

    Returns (Pi_a_new, Pi_b_new, p_a_new, p_b_new).
    """
    F_a, F_b, tau_a, tau_b = _forces_torques(
        x_a, x_b, q_a, q_b, I_a_body, I_b_body, m_a, m_b,
    )
    return (
        Pi_a + s * tau_a,
        Pi_b + s * tau_b,
        p_a  + s * F_a,
        p_b  + s * F_b,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Main integrator step
# ─────────────────────────────────────────────────────────────────────────────

def se3_splitting_step(
    state: np.ndarray,
    dt: float,
    m_a: float, I_a_body: np.ndarray,
    m_b: float, I_b_body: np.ndarray,
) -> np.ndarray:
    """One 4th-order Yoshida splitting step on SE(3) × SE(3).

    Sequence of 7 substeps:
      drift(c[0]·dt) → kick(d[0]·dt) → drift(c[1]·dt) → kick(d[1]·dt)
      → drift(c[2]·dt) → kick(d[2]·dt) → drift(c[3]·dt)

    Drop-in replacement for se3_variational_step with identical call signature.
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

    I_inv_a = np.linalg.inv(I_a_body)
    I_inv_b = np.linalg.inv(I_b_body)

    # Convert to canonical momenta
    Pi_a = I_a_body @ body_a.omega
    Pi_b = I_b_body @ body_b.omega
    p_a  = m_a * body_a.v
    p_b  = m_b * body_b.v
    q_a, q_b = body_a.q, body_b.q
    x_a, x_b = body_a.x, body_b.x

    # Yoshida 4-stage composition: drift(c) → kick(d) → ... → drift(c)
    for i in range(4):
        c = YOSHIDA_C[i]
        Pi_a, q_a, x_a = _drift_substep(Pi_a, q_a, x_a, p_a, I_inv_a, m_a, c * dt)
        Pi_b, q_b, x_b = _drift_substep(Pi_b, q_b, x_b, p_b, I_inv_b, m_b, c * dt)
        if i < 3:
            d = YOSHIDA_D[i]
            Pi_a, Pi_b, p_a, p_b = _kick_substep(
                Pi_a, Pi_b, p_a, p_b,
                x_a, x_b, q_a, q_b,
                I_a_body, I_b_body, m_a, m_b, d * dt,
            )

    # Convert back to (omega, v) and pack
    omega_a_new = I_inv_a @ Pi_a
    omega_b_new = I_inv_b @ Pi_b
    v_a_new = p_a / m_a
    v_b_new = p_b / m_b

    body_a_new = BodyState(x_a, v_a_new, q_a, omega_a_new)
    body_b_new = BodyState(x_b, v_b_new, q_b, omega_b_new)
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
    """Integrate n_steps with the Yoshida splitting and return trajectory.

    Returns a list of state arrays (including the initial state at index 0).
    """
    states = [state0]
    state = state0
    for i in range(n_steps):
        state = se3_splitting_step(state, dt, m_a, I_a_body, m_b, I_b_body)
        if (i + 1) % record_every == 0:
            states.append(state)
    return states
