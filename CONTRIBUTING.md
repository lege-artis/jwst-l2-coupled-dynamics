# Contributing to jwst-l2-coupled-dynamics

Thank you for your interest in contributing. This project is a numerical
reference implementation of two-rigid-body gravitationally-coupled dynamics
in the spirit of the **lege-artis** organisation — first-principles
derivations, primary-source citations, layered test pyramids, cross-backend
bit-identity where applicable.

Contributions in Python (reference layer), Fortran (Phase 1 canonical
backend), Pascal / C++ / Rust (planned Phase 2+ backends), or canonical-tier
mathematical documentation are welcome.

---

## Ground Rules

- All contributions must pass the project's **acceptance gates** before
  merging (see [Acceptance gates](#acceptance-gates))
- Code style must match the conventions already present in each backend
  directory (Python: PEP 8 + numpy docstrings; Fortran: free-form with
  `-std=f2018`, ASCII-only source per KB-039)
- One feature or fix per pull request — keep diffs reviewable
- Commit messages use the imperative: `Add Lagrange-top libration testbed`,
  `Fix gravity-gradient sign drift`
- Wikipedia is not a primary source. Citations to peer-reviewed textbook
  literature (Goldstein 2002, Marsden & Ratiu 1999, Arnold 1989,
  Hughes 2004, Thorne & Blandford 2017, Murray & Dermott 1999, Podolský 2018,
  etc.) are required for any new equation or physical claim.

---

## Acceptance gates

Every backend change must produce numerical outputs that match the Python
reference at element-wise tolerance per Higham 2002 §4.2 Option G
(ULP × sqrt(N)). The cross-test driver is at
`backends/fortran/tools/cross_test_python.sh` (Phase 1 Fortran reference);
the same pattern extends to additional backends.

The acceptance criteria are summarised in the project's
`docs/canonical/en/00-introduction.md` §4 (lege-artis doctrine §4.4
C1-C4 mapped onto this domain):

| Criterion | This project's interpretation |
|-----------|------------------------------|
| **C1** Textbook-verbatim test coverage | ≥3 analytical-limit configurations per canonical equation (Chapter 06 testbeds) |
| **C2** Cross-backend bit-identity | Fortran ≡ Python at element-wise ULP × sqrt(N) tolerance on shared fixtures |
| **C3** Per-routine textbook citation in source headers | grep-audited at PR review time |
| **C4** Empirical step-size scaling (substituted for FFT complexity) | RK4 energy drift scales as O(dt^4) on torque-free configurations |

Pull requests touching numerical code MUST attach test output showing C1-C2
gates remain green.

---

## Development Setup

### Python reference (always-available)

```bash
git clone https://github.com/lege-artis/jwst-l2-coupled-dynamics.git
cd jwst-l2-coupled-dynamics
python3 -m pip install numpy pytest                # only deps
python3 run_first_example.py                       # ~3-4 s, writes NDJSON
python3 -m pytest tests/                           # 43 tests, ~5 s wall
```

### Fortran reference backend (Phase 1)

Requires **gfortran** (`-std=f2018` support; gfortran 13+ tested) and **GNU
Make**. Cross-test additionally requires **Python 3 with NumPy** (the
fixture oracle).

```bash
cd backends/fortran
make clean && make all          # builds modules + runs all unit tests
make cross-test                 # runs C2 cross-test against Python fixtures
```

The Fortran backend uses ASCII-only source (per KB-039 lessons from
lege-artis/fourier), OS-aware Makefile (per KB-038), and an explicit
equation-to-code mapping in each module header.

### Pascal / C++ / Rust backends

Phase 2 work. Not in v0.1.0.

---

## Mathematical contributions

Canonical-tier chapters at `docs/canonical/en/` follow the standard set
in §2 of `docs/canonical/en/00-introduction.md`. Specifically:

- Equations have either a derivation in-chapter or a citation to a primary
  source where the derivation appears
- Coordinate-free formalism is preferred (configuration manifold SE(3) ×
  SE(3); Lie-group Euler-Poincaré reduction); specialisations to component
  form are flagged as such
- Two-in-parallel methodology where applicable (Approach A in Euler-angle
  chart + Approach B via Lie-group reduction, with cross-validation; see
  Ch 04 for the canonical example)
- Sign conventions are documented at the boxed-formula level (the
  OQ-FORT-1 sign-drift in v0.1.x development is the worked example —
  see Ch 03 §3.7.2)
- Open questions surface as `OQ-N.M` flags at chapter end and ride forward
  to v0.x.y resolution

---

## Bug Reports

Open an issue with:
- The exact `python3 run_first_example.py` or `make cross-test` output
  reproducing the issue
- Expected vs observed numerical values (or error output)
- Platform: OS, Python version, gfortran version (if applicable)
- Citation if the report concerns a discrepancy with a primary source

---

## Pull Request Checklist

- [ ] C1 gate: relevant analytical-limit testbeds still pass
- [ ] C2 gate: cross-backend bit-identity still holds (if touching numerical kernels)
- [ ] C3 gate: source-header citations updated for any new equation or
      modified algorithm
- [ ] Tests pass on `pytest tests/`
- [ ] Fortran backend (if touched) builds without warnings on at least one
      target platform
- [ ] PR description includes a one-sentence summary plus citation of the
      primary source motivating the change, if applicable

---

## Code of Conduct

This project follows the [Contributor Covenant](CODE_OF_CONDUCT.md).
All contributors are expected to uphold its standards.
