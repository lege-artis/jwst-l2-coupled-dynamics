# lege-artis/jwst-l2-coupled-dynamics — v0.1.0

**Released:** 2026-05-25
**Tag:** `v0.1.0`
**Status:** First public release. Python reference + Fortran Phase 1
canonical-reference backend; canonical-tier documentation Chapters 00..06
authored in English.

---

## What's in v0.1.0

This release establishes the **lege-artis** canonical-reference pattern for
two-rigid-body gravitationally-coupled dynamics — the third project in the
lege-artis numerical-reference family (after fourier and kh-sim), sharing
the same multi-backend, layered-test-pyramid, primary-source-anchored
discipline.

### Python reference layer (code, Apache 2.0)

The Python reference is the **oracle for all subsequent backends**. It
translates the canonical equations directly (not via any other source) and
serves as the C2 fixture provider.

| File | Purpose | Lines |
|------|---------|-------|
| `geometries.py` | Composite-body construction (cones/disks/cylinders/rods) + parallel-axis-theorem composition; `make_jwst_like()` prolate JWST-shaped reference body + `make_oblate_reference_body()` synthetic oblate for cross-stability testbed | ~430 |
| `dynamics.py` | Coupled 12-DOF equations: Newton (translation) + Euler (body-frame rotation) + quaternion kinematics. Newtonian gravity + gravity-gradient torque coupling. Classical RK4 integrator with explicit quaternion renormalisation. | ~290 |
| `run_first_example.py` | Worked example: tilted initial orientation + non-principal-axis spin → activates free precession + GG torque visibly. Writes NDJSON trajectory. | ~300 |
| `testcases.py` | 6-case L × direction sweep parameterisation | ~340 |
| `render_figures.py` | Diagnostic plots from NDJSON trajectory | ~190 |
| `tests/test_*.py` | 43 pytest tests across geometries, dynamics, conservation, regression | ~700 |
| `tests/fixtures/generate_fixtures.py` | Deterministic JSON oracle generator (3 fixtures: gravity_gradient_inputs, multipole_potential_inputs, symmetric_top_trajectory) | ~525 |

**Acceptance:** 43/43 pytest PASS in ~5 s wall. Worked example produces
`outputs/first_example_trajectory.ndjson` with conserved-quantity drift
`|dE/E₀| ≈ 6.2e-12` and `|dL/L₀| ≈ 2.1e-10` over 600 s at dt = 0.05 s
(machine-precision behaviour for classical RK4).

### Fortran reference backend (code, Apache 2.0)

Phase 1 backend covers what is closed by canonical-tier Chapters 01-03:
quaternions, geometries, gravity-gradient torque, multipole potential,
the orientation-response Jacobian, and a symmetric-top integrator suitable
for the closed-form testbed. Phase 2 (Lagrangian assembly, Euler-Poincaré
reduction, libration testbed, L2 frame) is the natural v0.2 scope.

| Module | Lines | Status |
|--------|-------|--------|
| `src/jwst_l2_constants.f90` | 45 | OQ-1 frozen physical constants (G_NEWTON to CODATA 2018) |
| `src/jwst_l2_quaternion.f90` | 125 | Line-for-line port of `dynamics.py:41-68` |
| `src/jwst_l2_geometries.f90` | 310 | All primitives + composite helpers, including `make_oblate_reference_body` |
| `src/jwst_l2_potential.f90` | (auth Phase A) | Mutual gravitational potential, multipole expansion |
| `src/jwst_l2_torque.f90` | (auth Phase A) | Gravity-gradient torque + cross-product |
| `src/jwst_l2_jacobian.f90` | (auth Phase A) | Orientation-response Jacobian (OQ-FORT-1 sign-fix landed) |
| `src/jwst_l2_diagnostics.f90` | (auth Phase A) | Energy + angular momentum conserved quantities |
| `src/jwst_l2_dynamics.f90` | (auth Phase A) | Symmetric-top integrator |

**13 unit tests** (test_quaternion, test_geometries with T1..T8 including
the oblate-body verification from OQ-FORT-6 Path C, test_potential,
test_torque, test_jacobian, test_diagnostics, test_dynamics,
test_anisotropy_sweep).

**3 C2 cross-tests** (Higham 2002 §4.2 Option G, tolerance = max(|ref|,
floor) × ULP × sqrt(N)):

| Test | Cases | Worst residual on v0.1.0 |
|------|-------|--------------------------|
| D4 gravity_gradient | 7 | 3.153e-16 (≈ 1 ULP) |
| D5 multipole_potential | 3 | 0.0 (exact) |
| D8 symmetric_top trajectory | 201 sampled states from 12000-step RK4 | 3.110e-11 at step 12000 |

All three C2 gates PASS. Tag landed as `jwst-l2-fortran-v0.1.0-phase-1-green`.

**Build chain (KB-038 + KB-039 from lege-artis/fourier):** OS-aware Makefile
(Windows cmd.exe + POSIX sh + macOS), ASCII-only source (gfortran
`-std=f2018 -pedantic` byte-counts UTF-8 against the 132-column limit),
case-aware identifier discipline.

### Canonical-tier documentation (CC-BY-SA-4.0)

| Chapter | Title | Lines | Status |
|---------|-------|-------|--------|
| 00 | Introduction | ~470 | Authored |
| 01 | Configuration manifold and variational principle | (scaffold + body) | Authored |
| 02 | Kinetic energy on SE(3) × SE(3) in intrinsic form | ~430 | Authored (with §2.5.1 d**I**/dt sensitivity, §2.6.1 CoM symplectic, §2.7.1 left-invariance) |
| 03 | Mutual gravitational potential and multipole expansion | ~640 | Authored (with §3.7.2 OQ-FORT-1 sign-fix and oblate/prolate stability contrast) |
| 04 | Euler-Lagrange equations and Euler-Poincaré reduction | ~440 | Authored (Approach A in Euler-angle chart + Approach B via Lie-group reduction, cross-validated) |
| 05 | Body-frame Euler equations as a specialisation | ~330 | Authored (principal-axis specialisation + symmetry-class worked paragraphs + libration cross-check) |
| 06 | Analytical limits and validation gates | ~325 | Authored (Kepler-limit + free symmetric top + GG libration + Dzhanibekov + decoupling testbeds + C2 procedural gate; with closed-form Lagrange-top (§6.4.3) and Jacobi-elliptic Dzhanibekov (§6.5.3) extensions) |

**Doctrine §4.4 mapping:**

| Criterion | This release |
|-----------|--------------|
| **C1** Textbook-verbatim coverage | ≥3 analytical-limit testbeds per canonical equation; 5 testbeds in Ch 06 with closed-form oracles |
| **C2** Cross-implementation bit-identity | Fortran ≡ Python at ULP × sqrt(N), 3 cross-tests on shared fixtures, all PASS |
| **C3** Per-routine textbook citation | Fortran module headers cite Goldstein/Hughes/Marsden-Ratiu/Landau-Lifshitz at routine granularity |
| **C4** Empirical scaling audit (substituted for FFT complexity) | RK4 energy drift O(dt⁴) — testbed at Ch 06 §6.7 |

### Reference bibliography

- Goldstein, Poole & Safko (2002) *Classical Mechanics* 3rd ed
- Arnold (1989) *Mathematical Methods of Classical Mechanics* 2nd ed
- Marsden & Ratiu (1999) *Introduction to Mechanics and Symmetry* 2nd ed
- Landau & Lifshitz (1976) *Mechanics* Vol I 3rd ed (Jacobi-elliptic
  asymmetric-top parametrisation §37)
- Thorne & Blandford (2017) *Modern Classical Physics* (rigid-body
  geometric framing + multipole expansion)
- Hughes (2004) *Spacecraft Attitude Dynamics* (gravity-gradient torque
  literature anchor)
- Murray & Dermott (1999) *Solar System Dynamics* (multipole expansion in
  the celestial-mechanics tradition)
- Jackson (1999) *Classical Electrodynamics* 3rd ed (multipole formalism)
- Podolský (2018) *Teoretická mechanika ve třech knihách* (Czech-language
  canonical anchor; Kapitoly 6-8 — transport equation, Euler kinematic/
  dynamic equations, Lagrangian formalism, free/damped/Lagrange symmetric tops)
- Higham (2002) *Accuracy and Stability of Numerical Algorithms* 2nd ed
  (Option G tolerance model for the C2 gate)
- Hairer, Lubich, Wanner (2006) *Geometric Numerical Integration* (forward
  reference for v0.2 Lie-group variational integrator)

### License stack

- Code: Apache License 2.0
- Documentation: CC-BY-SA-4.0
- Names: TRADEMARK.md declares MIM2000 / Improwave / Petr Yamyang protected

### Build-chain lessons applied (carried from lege-artis/fourier)

- **KB-037** — Fortran case-insensitivity (lowercase loop indices, `nlen`
  for size)
- **KB-038** — OS-aware Makefile (`ifeq ($(OS),Windows_NT)`)
- **KB-039** — ASCII-only Fortran source (gfortran `-pedantic` UTF-8 byte
  counting)
- **KB-043** — boxed-formula intra-doc inconsistency pattern caught by
  C2 cross-test discipline (the OQ-FORT-1 sign-fix is the worked example
  for v0.1.0)
- **KB-044** — bilingual canonical-tier CS/EN anchor convention
- **KB-045** — asymmetric rigid body conserved quantity E = ½ω·(Iω), not |ω|²
- **KB-046** — gravity-gradient instability measurement requires
  co-rotating spin + in-plane perturbation + last-quarter regression window
- **KB-047** — bash heredoc env export must precede the heredoc line
- **KB-048** — unit quaternion Higham tolerance uses scale_floor = 1.0

## Implementation discipline

- **First-principles translation:** Python and Fortran kernels both
  translate the canonical-tier equations directly. Neither is the
  source-of-truth for the other; `docs/canonical/en/*.md` is.
- **No external math libraries beyond NumPy:** the Python reference uses
  only NumPy primitives. Fortran uses only intrinsics. No external linear-
  algebra, ODE, or special-function libraries.
- **No optimisation in v0.x reference builds:** Fortran builds at `-O0
  -ffp-contract=off` (matches fourier v0.1.x convention). Optimised build
  paths are queued for v0.5+.
- **ASCII-only source.** Banner blocks use `=` × 57.
- **Equation-to-code mapping in module headers.** Every Fortran source
  file opens with a doc-comment block showing how each math operation
  maps to specific code lines.

## What's NOT yet in v0.1.0 (forward-looking roadmap)

- **Engineer-tier documentation (EN)** — quick-start, deploy/use, debug
  recipes. Queued v0.1.x.
- **CS canonical-tier translation** — Chapters 00-02 brief ready at
  `_handoffs/cs-canonical-tier-v0.1/SONNET-BRIEF.md` in the upstream
  workspace. Queued v0.1.1.
- **L2 restricted three-body extension** — rotating Earth-Sun frame,
  pseudo-potential, linearised stability at L2. Queued v0.2.
- **Symplectic Lie-group variational integrator** — Hairer-Lubich-Wanner
  2006 framework, alternative to RK4 with provable conservation. Queued v0.2.
- **Pascal / C++ / Rust backends** — mirroring lege-artis/fourier +
  lege-artis/kh-sim multi-backend layout. Queued v0.2.x.
- **Lyapunov-spectrum methods near Dzhanibekov regime** — quantify chaotic
  sensitivity in the unstable asymmetric-top dynamics. Queued v0.2.
- **Shad-tier engineer-narrated documentation band** — extending the
  existing `SHAD-NARRATIVE.md` to a full chapter sequence in the
  lege-artis/fourier band-structure pattern. Queued v0.3.
- **Visualisation tooling** — matplotlib → Plotly → Three.js per the
  4-tier visualisation plan. Queued v0.3.

## Citation

```
Yamyang, P. (2026). lege-artis/jwst-l2-coupled-dynamics v0.1.0 —
Python + Fortran canonical reference implementations of two-rigid-body
gravitationally-coupled dynamics. https://github.com/lege-artis/jwst-l2-coupled-dynamics
(tag v0.1.0)
```

## Acknowledgements

The discipline carried from lege-artis/fourier — first-principles equation
translation, validation against independent oracles, layered test pyramid,
ASCII-only source, OS-aware build chain, equation-to-code mapping in
headers — is what makes the second-language port land cleanly on first
attempt. Authoring the Fortran backend from `docs/canonical/en/*.md`
(rather than from `dynamics.py`) is what makes the cross-language match
meaningful: two independent renderings of the same math, both validated
against the same oracle, both producing numerically identical results.
That's not a coincidence; that's the canonical-reference pattern working
as designed.

The bilingual CS/EN canonical-tier anchor structure (KB-044) was developed
in this project's authoring round 3 and is now available as a reusable
pattern for any future bilingual lege-artis project.
