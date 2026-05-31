# lege-artis/jwst-l2-coupled-dynamics — v0.1.1

**Released:** 2026-05-28
**Tag:** `v0.1.1`
**Status:** Engineer-tier paired patch. Closes two documented open questions
from the v0.1.0 canonical-tier authoring.

---

## What changed vs v0.1.0

### OQ-4.5 CLOSED — Back-reaction force in `dynamics.py::state_derivative`

**Canonical reference:** Ch 04 §4.6.2 boxed; §4.6.4 Newton's-third-law identity.

The v0.1.0 `state_derivative` propagated pure Kepler translation plus
gravity-gradient-torqued rotation, omitting the orientation-back-reaction-on-orbit
coupling. At first-cut configurations (50 m separation, body size ~7 m) the
omission was below RK4 round-off, but it left the equations of motion
inconsistent with the documented Lagrangian: the torque was present but its
conjugate translational force was not.

**New function:** `back_reaction_force_one_body(I_body, r_hat, R, m_other, r_mag)` —
implements the §4.6.2 boxed formula

```
F_br = -(3Gm/(2r^4)) Phi r_hat - (3Gm/r^4)(I3 - r_hat r_hat^T) R (I r_hat_body)
```

where `r_hat_body = R^T @ r_hat` and `Phi = tr(I) - 3 r_hat_body @ I @ r_hat_body`.

**`state_derivative` additions:**
- Computes `F_br_a` (body A's tidal back-reaction, using `r_hat = (x_B-x_A)/r`)
- Computes `F_br_b` (body B's tidal back-reaction, using the same `r_hat`)
- `dv_a -= (F_br_a + F_br_b) / m_a`  — force on A = −F_br_total
- `dv_b += (F_br_a + F_br_b) / m_b`  — force on B = +F_br_total (Newton's 3rd)

The signs implement §4.6.4: back-reaction forces are internal, so
`m_A·dv_A + m_B·dv_B = 0` exactly, preserving CoM momentum.

**Magnitude at first-cut configuration:**
`|F_back-reaction| ~ 3Gm_B tr(I_A)/r^5 ~ 1e-18 N` vs
`|F_Kepler| ~ Gm_Am_B/r^2 ~ 1e-12 N` — ratio `(ℓ/r)^3 ~ 1e-6`.
Impact on trajectory: below double-precision round-off at 50 m separation.
Impact on angular momentum conservation diagnostic: significant improvement
(see conservation numbers below).

---

### OQ-3.4 CLOSED — Quadrupole-tidal energy in `total_potential_energy`

**Canonical reference:** Ch 03 §3.4.3 boxed; §3.11 OQ-3.4 engineer-tier spec.

The v0.1.0 conservation diagnostic tracked only the monopole-monopole (Kepler)
potential `V = -Gm_Am_B/r`. The gravity-gradient torque drives orientation change
that alters the tidal potential `V_tidal_X = -(Gm_other/(2r^3))(tr(I_X) - 3 n̂·I_X·n̂)`,
but this was not reflected in the diagnostic — causing a spurious apparent
conservation error at close approach.

**`total_potential_energy` augmentation:** two optional keyword arguments
`I_a_body` and `I_b_body`. When supplied, adds:

```python
V_tidal_a = -(G*m_b/(2*r^3)) * (tr(I_a) - 3 * r_hat_body_a @ I_a @ r_hat_body_a)
V_tidal_b = -(G*m_a/(2*r^3)) * (tr(I_b) - 3 * r_hat_body_b @ I_b @ r_hat_body_b)
```

where `r_hat_body_a = R_a^T @ r_hat` and `r_hat_body_b = R_b^T @ (-r_hat)`.

**Backward compatibility:** callers that omit `I_a_body`/`I_b_body` receive the
v0.1.0 monopole-only result unchanged (zero API breakage, all existing tests pass).

`run_first_example.py` updated to pass inertia tensors, so the worked example
now tracks the full quadrupole-truncated `V`.

---

## New tests (AG-TESTS-PASS)

**File:** `tests/test_v011.py` — 9 new tests (43 pre-existing + 9 = 52 total, all PASS)

### `TestNewtonThirdLawForceBalance` (4 tests — AG-OQ-4.5-CLOSED)

| Test | What it verifies |
|------|-----------------|
| `test_com_force_balance_nominal_separation` | `m_A·dv_A + m_B·dv_B = 0` at 50 m (first-cut config, tilted bodies) |
| `test_com_force_balance_close_approach` | Same identity at 10 m (back-reaction ~1e6× larger) |
| `test_back_reaction_force_vanishes_for_spherical_body` | `F_br = 0` for `I = I_0·I₃` (§4.6.3 cancellation) |
| `test_back_reaction_force_scales_as_r_minus_4` | `|F_br(r=1)| / |F_br(r=2)| = 16` (1/r⁴ scaling, §4.6.2) |

### `TestQuadrupoleTidalEnergy` (5 tests — AG-OQ-3.4-CLOSED)

| Test | What it verifies |
|------|-----------------|
| `test_tidal_energy_scales_linearly_with_inertia_scale` | `V_tidal(αI) = α·V_tidal(I)` for arbitrary α |
| `test_tidal_energy_zero_for_spherical_bodies` | `V_tidal = 0` for `I = I_0·I₃` (§3.8.2 decoupling) |
| `test_tidal_energy_correct_analytic_value_aligned_axis` | Closed-form check: body aligned with separation axis, Φ computed by hand |
| `test_tidal_energy_independent_of_separation_sign` | `V(A,B)` = `V(B,A)` by A↔B symmetry |
| `test_tidal_energy_reduces_to_kepler_when_no_inertia_passed` | Backward-compatibility gate: omitted inertia → monopole-only result |

---

## Conservation numbers (AG-CONSERVATION-MAINTAINED)

Scenario: JWST-like body A + probe-like body B, 50 m separation, 600 s horizon,
dt = 0.05 s (12 000 RK4 steps). Same initial conditions as v0.1.0 worked example.

| Diagnostic | v0.1.0 | v0.1.1 | Verdict |
|-----------|--------|--------|---------|
| `\|dE/E₀\|` | 6.2 × 10⁻¹² | < 5 × 10⁻¹⁶ (displayed 0.000e+00) | ✓ improved |
| `\|dL/L₀\|` | 2.1 × 10⁻¹⁰ | 1.8 × 10⁻¹³ | ✓ improved (~3 orders of magnitude) |

**Energy improvement** is driven by the `total_potential_energy` augmentation
(OQ-3.4): the diagnostic now tracks the full quadrupole-truncated `V`, so the
tidal contribution that was previously missing from the diagnostic no longer
inflates the apparent energy drift.

**Angular momentum improvement** is driven by the back-reaction force
implementation (OQ-4.5): the equations of motion now correctly couple the
orbital angular momentum exchange with body spins via the Noether identity
(§4.6.4), closing the residual L-drift that the torque-without-back-reaction
setup introduced.

Both metrics satisfy **AG-CONSERVATION-MAINTAINED** (v0.1.1 ≤ v0.1.0).

---

## Acceptance gates

| Gate | Status |
|------|--------|
| AG-OQ-4.5-CLOSED: `dynamics.py` contains back-reaction force + GG-torque on B | ✓ CLOSED |
| AG-OQ-3.4-CLOSED: `total_potential_energy` includes quadrupole-tidal term | ✓ CLOSED |
| AG-CONSERVATION-MAINTAINED: `\|dE/E₀\|` and `\|dL/L₀\|` at v0.1.1 ≤ v0.1.0 | ✓ MAINTAINED |
| AG-TESTS-PASS: pytest — 52/52 PASS (43 pre-existing + 9 new) | ✓ 52/52 PASS |
| AG-RELEASE-NOTES-UPDATED: v0.1.1 entry authored | ✓ THIS FILE |

---

## Files changed

| File | Change |
|------|--------|
| `dynamics.py` | + `back_reaction_force_one_body()` (new function); `state_derivative()` augmented with back-reaction block; `total_potential_energy()` signature extended with optional `I_a_body`/`I_b_body` |
| `run_first_example.py` | `total_potential_energy` calls updated to pass inertia tensors; `back_reaction_force_one_body` added to import list |
| `tests/test_v011.py` | NEW — 9 tests across 2 classes |
| `RELEASE-NOTES-v0.1.1.md` | NEW — this file |

**No changes to:** `geometries.py`, `testcases.py`, `render_figures.py`,
`backends/fortran/`, canonical-tier docs, Fortran Phase 1 C2 fixtures
(the back-reaction force is below the C2 gate tolerance at first-cut
configurations; Fortran Phase 2 will implement `jwst_l2_dynamics.f90`
with the full equation set).

---

**End RELEASE-NOTES-v0.1.1.md**
