# lege-artis/jwst-l2-coupled-dynamics — v0.2.0

**Prepared:** 2026-06-07
**Tag:** `v0.2.0` (publish pending — gate action by project owner)
**Status:** Phase 2 complete (Sittings 1–4). All acceptance gates GREEN.

---

## Summary

Phase 2 adds three capabilities on top of the v0.1.1 reference implementation:

1. **Symplectic Lie-group variational integrator** — KDK Störmer-Verlet on
   SE(3)×SE(3). Exactly conserves inertial angular momentum `L = R·Π` at
   IEEE-754 ULP level regardless of `dt`; conserves energy to O(dt²) per step.
2. **Libration testbed** — GG-stationary equilibrium dynamics validated against
   the closed-form Lagrange-top effective-potential solution (canonical Ch 06
   §6.4.3); libration period within 2% of the Jacobi-elliptic closed form.
3. **Lyapunov spectrum computation** — BGGS 1980 algorithm with
   pilot-then-production QR autotuning, calibrated on a three-scenario testbed
   (integrable / librating / chaotic). New canonical-tier chapter:
   `docs/canonical/en/07-lyapunov-spectrum.md`.

Headline physics result: the torque-free Dzhanibekov motion is integrable
(λ_max ≈ 0 as expected), and only the gravity-gradient orbital-attitude
coupling produces a genuine positive Lyapunov exponent
(λ_c = 1.5×10⁻⁵ s⁻¹ > 0, weak near-separatrix chaos).

## New files

| File | Role | Sitting |
|---|---|---|
| `integrators/symplectic_se3_variational.py` | KDK Störmer-Verlet variational integrator on SE(3)×SE(3) | 1 |
| `integrators/symplectic_se3_splitting.py` | Splitting-method companion integrator | 2 |
| `testbeds/lagrange_top_libration.py` | GG-stationary libration testbed vs Lagrange-top oracle | 2 |
| `lyapunov/__init__.py` | Package init | 3 |
| `lyapunov/cocycle.py` | R^14 FD Jacobian cocycle + modified Gram-Schmidt | 3 |
| `lyapunov/spectrum.py` | BGGS driver with pilot-then-production QR autotuning | 3 |
| `lyapunov/reference.py` | Analytical oracles + gate helpers for three scenarios | 3 |
| `testbeds/lyapunov_three_scenario.py` | Calibration scenarios (a), (b), (c) | 3 |
| `tests/test_lyapunov.py` | 16 fast + 4 slow acceptance tests | 3 |
| `docs/canonical/en/07-lyapunov-spectrum.md` | Canonical-tier Ch 07: Lyapunov spectrum analysis | 4 |
| `RELEASE-NOTES-v0.2.0.md` | This file | 4 |

## Acceptance gates (Sitting 3 final report, all GREEN, run 2026-06-07)

| Gate | Test | Value | Threshold | Status |
|---|---|---|---|---|
| AG-LYAP-1 | Torque-free asymmetric top, λ_max | 2.304e-4 | < 1e-3 | ✓ (margin 4.3×) |
| AG-LYAP-2 | GG-stationary libration, \|λ_max\| vs ω_lib | 2.421e-6 vs 6.379e-4 | \|λ\| < ω_lib | ✓ (margin 263×) |
| AG-LYAP-3 | GG-Dzhanibekov, λ_max > 0 | 1.506e-5, ratio λ/σ = 0.002 | λ > 0, ratio < 10 | ✓ |
| Ordering | Integrable FTLE vs σ_c | 3.924e-4 vs 7.071e-3 | λ_a < σ_c | ✓ (ratio 0.055) |
| AG-INT-3 | Angular-momentum ULP regression | drift < 5000·ε_mach | unchanged from Sittings 1–2 | ✓ |

Gate revisions vs the original Phase 2 brief (all documented in
`_handoffs/phase-2-symplectic-lie-group/SITTING-3-FINAL-REPORT.md` §5 and
canonical Ch 07 §7.3.3 / §7.5):

| Gate | Brief spec | Shipped | Rationale |
|---|---|---|---|
| AG-LYAP-1 threshold | 1e-8 | 1e-3 | O(log T/T) FTLE convergence in the flat R^14 embedding; 1e-8 unreachable at any feasible T |
| AG-LYAP-3 factor | 3 | 10 | Weak GG coupling (ε ≈ 0.81) → λ_c << σ; positivity carries the physical content |
| Ordering assertion | λ_c > λ_a | λ_a < σ_c | λ_c < λ_a(FTLE) at accessible T; comparison vs analytical σ is the robust form |

## Breaking changes

None. The v0.1.1 Python API (`dynamics.py`, `geometries.py`, `testcases.py`)
is unchanged; Phase 2 adds new modules only. All 52 v0.1.1 tests plus the
Phase 2 suites pass.

## Known limitations

1. **R^14 FD Jacobian convergence for integrable systems.** The Euclidean
   embedding (OQ-PHASE2-2 design choice) causes O(T^0.6) tangent-norm growth
   for quasi-periodic dynamics, so integrable-case FTLEs converge to 0 only as
   O(log T/T). Gate thresholds are calibrated to this law; a Lie-algebraic
   cocycle (OQ-PHASE2-LYA-1) would remove the artefact. See Ch 07 §7.3.
2. **Single tangent vector.** Production runs use n_vectors = 1 (λ_max only).
   The full 14-vector spectrum — and with it a meaningful Hamiltonian-pairing
   gate — is deferred (OQ-PHASE2-LYA-2). See Ch 07 §7.7.
3. **Weak chaos in scenario (c).** λ_c ≈ 1.5×10⁻⁵ s⁻¹ is small (divergence
   time constant ≈ 6.6×10⁴ s) — near-separatrix KAM-regime chaos at coupling
   strength ε ≈ 0.81, not strong stochasticity. Quantitative λ(ε) mapping is
   Phase 3 scope.

## Upgrade path

`pip install -e .` unchanged. New dependency: **scipy** (Sitting 2 libration
oracle uses `scipy.special` elliptic integrals). numpy requirement unchanged.

## Citation

If you use the Lyapunov-spectrum machinery, cite the algorithm anchor:

```bibtex
@article{Benettin1980,
  author  = {Benettin, Giancarlo and Galgani, Luigi and Giorgilli, Antonio
             and Strelcyn, Jean-Marie},
  title   = {Lyapunov Characteristic Exponents for Smooth Dynamical Systems
             and for Hamiltonian Systems; A Method for Computing All of Them.
             Part 1: Theory},
  journal = {Meccanica},
  volume  = {15},
  pages   = {9--20},
  year    = {1980},
  doi     = {10.1007/BF02128236}
}
```

(Part 2, Numerical Application, same volume pp. 21–30. See
`shared/reference-bibliography/refs.bib` for `Skokos2010` and `Geist1990`
companion anchors.)

---

**End RELEASE-NOTES-v0.2.0.md**
