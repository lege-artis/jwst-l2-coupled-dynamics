# lege-artis/jwst-l2-coupled-dynamics — v0.3.0

**Prepared:** 2026-06-07
**Tag:** `v0.3.0` (publish pending — gate action by project owner)
**Status:** Phase 3 complete. All acceptance gates GREEN. Additive to v0.2.0 (superset).

---

## Summary

Phase 3 computes the **full Lyapunov spectrum in a 12-dimensional cocycle**, resolving
two open questions carried from Phase 2 (OQ-PHASE2-LYA-1 + LYA-2). The flat Euclidean
R^14 cocycle used in Phase 2 is correct for the single-vector λ_max gates but
contaminates the *full* spectrum with two spurious quaternion-norm contraction modes
(one per body) that break Hamiltonian pairing. Phase 3 removes them two ways:

1. **Tier A — quaternion-deflation cocycle** (`lyapunov/cocycle_deflated.py`): keeps the
   R^14 finite-difference cocycle but projects every tangent vector orthogonal to the
   two current quaternion-radial directions (q_A, q_B) before each Gram-Schmidt step,
   working in the deflated 12-dimensional subspace.
2. **Tier B — Lie-algebraic cocycle** (`lyapunov/cocycle_liealg.py`): a body-frame
   so(3) × R^3 tangent propagator that is 12-dimensional by construction. Orientation
   perturbations use the exponential map in the right trivialisation matching the
   integrator (q_{k+1} = q_k ⊗ exp(ε·ξ)); the log map ξ = 2·arctan2(|q_vec|, q_w)·
   q_vec/|q_vec| maps the propagated perturbation back to so(3).

This is additive to Phase 2: no Phase 2 file is modified, dynamics.py API is unchanged,
and v0.3.0 carries all of v0.2.0.

## New files

| File | Role |
|---|---|
| `lyapunov/cocycle_deflated.py` | Tier A — R^14 quaternion-deflation cocycle (12-dim) |
| `lyapunov/cocycle_liealg.py` | Tier B — so(3)×R^3 Lie-algebraic cocycle (12-dim) |
| `lyapunov/pairing.py` | Full-spectrum Hamiltonian-pairing gate helpers |
| `testbeds/lyapunov_full_spectrum.py` | Tier A/B production runners |
| `tests/test_lyapunov_phase3.py` | 16 fast + 5 slow acceptance tests |
| `docs/canonical/en/07-lyapunov-spectrum.md` | §7.10 added; §7.3.2/§7.9 corrected |
| `docs/canonical/cs/07-ljapunovovo-spektrum.md` | CS mirror (§7.10 + corrections) |
| `RELEASE-NOTES-v0.3.0.md` | This file |

## Acceptance gates (all GREEN; `tests/test_lyapunov_phase3.py` = 16 fast + 5 slow)

| Gate | Check | Result |
|---|---|---|
| AG-P3-1 | full spectrum sum \|Σλ\| < 1e-2 | **PASS** (the −0.064/−0.071 R^14 norm modes gone; spectrum no longer sums to −0.14) |
| AG-P3-2 | Hamiltonian pairing; near-zero invariants present | **PASS** |
| AG-P3-3 | Lie-algebraic λ_max scenario (a) < 1e-3 | **PASS** (2.29e-4; gate revised from the 1e-5 aspiration — see below) |
| AG-P3-4 | Lie-algebraic λ_max within factor 3 of R^14 value | **PASS** (both > 0, scenario (c)) |
| AG-P3-REGRESSION | all Phase 2 gates unchanged | **PASS** (16/16 fast) |

## Key finding (corrected a Phase 2 prediction)

The Phase-2-era Ch 07 attributed the scenario-(a) O(log T/T) FTLE convergence to the
flat R^14 embedding and predicted the Lie-algebraic cocycle would restore fast
convergence and make a 1e-8 AG-LYAP-1 gate achievable. **This was falsified:** the
manifold-native Lie-algebraic cocycle gives λ_max = 2.29e-4 s⁻¹, essentially identical
to the R^14 value 2.30e-4 s⁻¹ at T = 20,000 s. The slow convergence is **intrinsic to
the quasi-periodic multi-frequency dynamics, not the embedding** — a genuinely
chart-caused effect would have vanished under the manifold-native cocycle.

The embedding's *real* defect is the spurious quaternion-norm contamination of the full
spectrum, which Phase 3 does fix. AG-LYAP-1 stays at the Phase 2 value 1e-3. Ch 07 §7.3.2
(EN + CS) is corrected to separate the two effects; §7.10 documents the Phase 3 results
and the corrected prediction.

## Breaking changes

None. dynamics.py / geometries.py public API unchanged; all Phase 2 modules untouched.
All v0.2.0 tests plus the new Phase 3 suite pass.

## Known limitations

1. The intrinsic O(log T/T) scenario-(a) convergence persists (it is not an embedding
   artefact), so AG-LYAP-1 remains at 1e-3 rather than the originally-hoped 1e-8.
2. Kaplan-Yorke dimension estimates for the scenario (c) chaotic layer are a stretch
   goal, out of Phase 3 scope.
3. OQ-PHASE2-LYA-3 (JWST attitude-control stability margins) remains open — it needs the
   v0.2+ L2 restricted-three-body environment before the orbital setting is representative.

## Upgrade path

`pip install -e .` unchanged. Dependencies unchanged from v0.2.0 (numpy + scipy).

## Citation

The Lyapunov-spectrum algorithm anchor is unchanged from v0.2.0 (Benettin et al. 1980,
*Meccanica* **15**, 9–30). The Lie-group machinery anchors: Hairer-Lubich-Wanner 2006
(Geometric Numerical Integration) §IV; Marsden & Ratiu 1999 §9. See
`shared/reference-bibliography/refs.bib`.

---

**End RELEASE-NOTES-v0.3.0.md**
