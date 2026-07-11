# lege-artis / jwst-l2-coupled-dynamics

[![Apache 2.0](https://img.shields.io/badge/code-Apache%202.0-blue.svg)](LICENSE)
[![CC BY-SA 4.0](https://img.shields.io/badge/docs-CC%20BY--SA%204.0-green.svg)](LICENSE-DOCS)
[![Citation](https://img.shields.io/badge/cite-CITATION.cff-orange.svg)](CITATION.cff)

Canonical reference implementation of two-rigid-body gravitationally-coupled
dynamics, anchored on a JWST-like primary and a probe-like secondary
satellite in deep space. Sibling project to
[lege-artis/fourier](https://github.com/lege-artis/fourier) and
lege-artis/kh-sim — same multi-backend, layered-test-pyramid pattern.

## Status

**v0.3.1 — doc-patch release** (2026-07): corrects the Ch 06 §6.5.3 boxed
Jacobi-elliptic parametrisation (the published v0.1.0–v0.3.0 boxed forms for
the amplitudes, characteristic frequency Ω\* and modulus k² were erroneous —
see the dated correction note in `docs/canonical/en/06-analytical-limits-and-validation-gates.md`
§6.5.3 and the new regression guard `tests/test_elliptic_oracle.py`, which
validates the corrected closed form against a tight-tolerance integration to
≲1e-12 in period and 8e-14 element-wise). Also refreshes this README, which
had lagged two releases.

Cumulative shipped state (v0.1.0 → v0.3.1):

- **Phase 2 (shipped inside v0.3.0)** — symplectic KDK Störmer-Verlet
  integrator on SE(3)×SE(3) (`integrators/`), Lagrange-top libration testbed
  with the (2/π)K(k) finite-amplitude oracle (`testbeds/`), BGGS
  Lyapunov-spectrum machinery with three calibrated scenarios (`lyapunov/`),
  canonical-tier Ch 07 (EN + CS).
- **Phase 3 (v0.3.0, 2026-07-04)** — Lie-algebraic so(3)×R³ cocycle +
  quaternion-deflation tier for the full 12-dimensional Lyapunov spectrum;
  Hamiltonian-pairing gate; Ch 07 §7.10 corrected-prediction narrative (the
  slow FTLE convergence is intrinsic to the quasi-periodic dynamics, not the
  R¹⁴ embedding).
- **Canonical-tier docs** — EN Ch 00–07 complete; **CS Ch 00–07 complete**
  (Podolský 2018 primary anchor); engineer-tier Ch 07 how-to (EN + CS).

The v0.1.0 foundation (2026-05-25):

- **Python reference layer** — runnable in ~3 s on a modern laptop. Coupled
  12-DOF dynamics (Newton for translation + Euler for rotation, body-frame),
  quaternion-based attitude, classical RK4, no external dependencies beyond
  NumPy. 43-test pytest suite. Demonstrates free precession at the Euler
  frequency and conservation of total energy and angular momentum at
  machine precision.
- **Fortran reference backend (Phase 1)** — green at tag
  `jwst-l2-fortran-v0.1.0-phase-1-green`. Builds with gfortran 13+
  (`-std=f2018`), ASCII-only source, OS-aware Makefile. 13 unit tests +
  3 C2 cross-tests against the Python first-cut at ULP × sqrt(N) tolerance
  (Higham 2002 §4.2 Option G). The C2 discipline catches sign-convention
  drift, off-by-one errors in cross-product order, and similar physical-
  truth-anchored bugs that documentation review alone misses.
- **Canonical-tier documentation (English)** — Chapters 00 through 06
  authored. Coverage:
  - Configuration manifold on SE(3) × SE(3) and the variational principle
  - Kinetic energy in intrinsic form (with d**I**/dt sensitivity bridge)
  - Mutual gravitational potential and multipole expansion (Cartesian +
    spherical-harmonic, with cross-validation at quadrupole order)
  - Euler-Lagrange equations and Euler-Poincaré reduction (Approach A in
    Euler-angle chart, Approach B via Lie-group reduction, cross-validated)
  - Body-frame Euler equations as a specialisation
  - Analytical limits and validation gates (5 testbeds, including
    closed-form Lagrange-top libration and Jacobi-elliptic Dzhanibekov
    parametrisation)
- **First-cut prototype trajectory** — `outputs/first_example_trajectory.ndjson`
  is a worked example demonstrating free precession at the predicted Euler
  frequency, with FFT-detected period matching theory to within one bin.

## Scope

| Domain element | v0.1 treatment |
|----------------|---------------|
| Two rigid bodies (extended, full inertia tensors) | Yes |
| Mutual gravitational potential | Yes, derived from first principles via multipole expansion |
| Rigid rotation about centre of mass | Yes, derived via Lagrange + Euler-Poincaré reduction |
| Relative orbital motion in centre-of-mass frame | Yes |
| Restricted three-body L2 dynamics (Sun-Earth-spacecraft) | Planned v0.4.x (CR3BP TG-8 oracles specified) |
| Relativistic 1PN/2PN corrections | Explicitly excluded (v/c ~ 10^-5 for JWST scales) |
| Solar radiation pressure, third-body Moon | Out of canonical-tier scope |
| Internal flexibility / vibration modes | Out of scope (rigid-body baseline) |

## Plan

* **v0.1.0** *(2026-05-25)* — Python reference + Fortran Phase 1 + canonical-tier Ch 00..06 (EN) — **shipped**
* **v0.1.1** *(2026-05-31)* — back-reaction force + quadrupole-tidal energy; conservation to ULP ceiling — **shipped**
* **v0.3.0** *(2026-07-04)* — symplectic Lie-group variational integrator +
  libration testbed + full 12-dim Lyapunov spectrum (Phases 2+3); CS canonical
  tier — **shipped** (supersedes the former v0.2.0 line)
* **v0.3.1** *(this release)* — Ch 06 §6.5.3 erratum + elliptic-oracle
  regression guard + README refresh — **shipped**
* **Next: v0.4.x** — L2 restricted three-body extension (CR3BP rotating frame;
  Master-Test-Set TG-8 oracles already specified); Master Test Set gate suite
  integration
* **Later** — Multi-backend ports: Pascal (NR 1986 anchor) + C++ + Rust,
  mirroring the lege-artis/fourier + lege-artis/kh-sim layout; Shad-tier
  engineer-narrated docs; visualisation tooling (matplotlib → Plotly →
  Three.js per tier plan)

## How to run

### Python reference (minimal)

```bash
git clone https://github.com/lege-artis/jwst-l2-coupled-dynamics.git
cd jwst-l2-coupled-dynamics
python3 -m pip install numpy pytest      # only deps
python3 run_first_example.py             # full 600 s coupled-dynamics simulation
python3 -m pytest tests/ -v              # 43 tests, ~5 s wall
```

Output: `outputs/first_example_trajectory.ndjson` (~1 MB; 601 NDJSON
snapshots). Each line is a JSON object with the full state of both
bodies — positions, velocities, quaternions, rotation matrices, body-frame
and inertial-frame angular velocities, inertia tensors in both frames,
plus conserved-quantity and torque-magnitude diagnostics.

### Fortran reference backend (Phase 1)

```bash
cd backends/fortran
make clean && make all      # 13 unit tests
make cross-test             # 3 C2 cross-tests vs Python fixtures
```

Requires gfortran 13+ and Python 3 with NumPy. The cross-test driver
compares Fortran outputs against `tests/fixtures/*.json` (deterministic
Python oracle) element-wise at Higham §4.2 Option G tolerance.

## What this demonstrates

1. **Composite inertia tensors with non-trivial principal moments**, built
   from cones, disks, cylinders, and rods via the parallel-axis theorem.
   Body A is JWST-shaped (sunshield + bus + boom + primary mirror,
   approximated as a prolate composite); body B is a probe (cylinder +
   dish antenna).

2. **Free precession at the Euler frequency.** Both bodies are given
   initial angular velocities with both a perpendicular component and a
   parallel (body-z) component. Euler's equation
   `I·dω/dt + ω × (I·ω) = τ` produces perpendicular-component oscillation
   at frequency `λ = (I_zz − I_xx) / I_xx · ω_z`. **Observed in FFT of
   the trajectory:** peak matches theory within one FFT bin.

3. **Conservation laws as numerical-method audit.** Without external
   bodies, total energy and total angular momentum must be conserved.
   At dt = 0.05 s over 600 s of integration:
   - `|dE/E₀| ≈ 6 × 10⁻¹²` (essentially machine precision)
   - `|dL/L₀| ≈ 2 × 10⁻¹⁰`
   This is the kind of validation that catches integrator bugs early.

4. **Gravity-gradient torque magnitude.** At 50 m separation, body A
   experiences a gravity-gradient torque from body B of order 10⁻⁹ N·m.
   For comparison, A's rotational KE at ω ≈ 0.08 rad/s is ~50 J. The
   GG torque modulates the free precession on a long timescale (~10⁶ s);
   the v0.1 worked example captures only the free-precession dominant
   mode.

5. **Cross-backend bit-identity (C2 gate).** Python and Fortran
   implementations are independently translated from
   `docs/canonical/en/*.md`, then compared element-wise on shared JSON
   fixtures. Worst-case observed residuals on the v0.1.0 green-tag run:
   - D4 gravity-gradient torque: 3.15e-16 (≈ 1 ULP, 7 cases)
   - D5 multipole potential: 0.0 (exact, 3 cases)
   - D8 symmetric-top trajectory: 3.11e-11 at step 12000 (201 sampled states)

## Documentation tiers

| Tier | Status | Audience | Location |
|------|--------|----------|----------|
| **Canonical (EN)** | Ch 00..07 complete (incl. §6.5.3 correction note 2026-07-11) | Mathematician, physicist, academic reviewer | `docs/canonical/en/` |
| **Canonical (CS)** | Ch 00..07 complete; Podolský 2018 primary anchor | Czech-speaking reader | `docs/canonical/cs/` |
| **Engineer (EN + CS)** | Ch 07 Lyapunov how-to shipped; remaining chapters queued | Practitioner who wants to use the code | `docs/engineer/` |
| **Shad-tier** | Single narrative shipped | Engineer-narrated walkthrough; oscilloscope-curious reader | `SHAD-NARRATIVE.md` |

## License

* **Code** (everything under `backends/`, `tests/`, Makefiles, `tools/`,
  `*.py`, GitHub Actions): Apache License 2.0 — see `LICENSE`
* **Documentation** (everything under `docs/`, `shared/`, `SHAD-NARRATIVE.md`,
  `RELEASE-NOTES-*.md`, this file): Creative Commons Attribution-ShareAlike
  4.0 International — see `LICENSE-DOCS`
* **Trademark notice**: `MIM2000`, `Improwave`, and the personal name
  `Petr Yamyang` are not licensed for derivative use — see `TRADEMARK.md`

## Central references

For analytical-mechanics formalism: Goldstein 2002, Arnold 1989, Marsden
& Ratiu 1999, Thorne & Blandford 2017, Landau-Lifshitz Vol I 1976.
For spacecraft attitude dynamics: Hughes 2004. For gravitational
potentials of extended bodies: Murray & Dermott 1999, Jackson 1999
(multipole formalism transferred from electrodynamics). For Czech-language
canonical anchor: Podolský 2018 *Teoretická mechanika ve třech knihách*,
Kapitoly 6-8. Full bibliography in `shared/reference-bibliography/refs.bib`.

## Citation

```
Yamyang, P. (2026). lege-artis/jwst-l2-coupled-dynamics v0.1.0 —
Python + Fortran canonical reference implementations of two-rigid-body
gravitationally-coupled dynamics. https://github.com/lege-artis/jwst-l2-coupled-dynamics
(tag v0.1.0)
```

## See also

* Sibling: [lege-artis/fourier](https://github.com/lege-artis/fourier) —
  DFT/FFT canonical reference; consumed by this project's tumble-spectrum
  diagnostic
* Sibling: lege-artis/kh-sim — Kelvin-Helmholtz 2D incompressible
  Navier-Stokes reference; same multi-backend pattern
* Doctrine: shared reference-implementation acceptance criteria
  (C1 textbook-verbatim tests, C2 cross-backend bit-identity, C3
  per-routine source citations, C4 empirical complexity audit)
