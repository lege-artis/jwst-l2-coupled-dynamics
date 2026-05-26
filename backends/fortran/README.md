# jwst-l2-coupled-dynamics — Fortran reference backend (Phase 1)

This directory holds the lege-artis Fortran canonical-reference implementation
of the two-rigid-body coupled-dynamics primitives. Phase 1 scope covers what
is closed by canonical-tier Chapters 01-03; Phase 2 (Lagrangian assembly,
Euler-Poincaré reduction, libration testbed, L2 frame) is queued for v0.2.

## Status

**Phase 1 green** at the v0.1.0 release tag.
13 unit tests + 3 C2 cross-tests PASS against Python fixtures at Higham §4.2
Option G tolerance. See `_audit/cross_test_summary.md` for the latest
cross-test report.

## Layout

```
backends/fortran/
├── README.md             — this file
├── Makefile              — OS-aware build + test driver (KB-038 pattern)
├── src/
│   ├── jwst_l2_constants.f90       — frozen physical constants
│   ├── jwst_l2_quaternion.f90      — quaternion algebra
│   ├── jwst_l2_geometries.f90      — primitives + composite-body helpers
│   ├── jwst_l2_potential.f90       — multipole potential
│   ├── jwst_l2_torque.f90          — gravity-gradient torque
│   ├── jwst_l2_jacobian.f90        — orientation-response Jacobian
│   ├── jwst_l2_diagnostics.f90     — energy + angular momentum
│   └── jwst_l2_dynamics.f90        — symmetric-top integrator
├── tests/
│   ├── test_quaternion.f90
│   ├── test_geometries.f90         — T1..T8 (includes oblate-body check)
│   ├── test_potential.f90
│   ├── test_torque.f90
│   ├── test_jacobian.f90           — T2a libration-eigenvalue regression
│   ├── test_diagnostics.f90
│   ├── test_dynamics.f90
│   └── test_anisotropy_sweep.f90
├── tools/
│   └── cross_test_python.sh        — Python-oracle cross-test runner (C2 gate)
└── _audit/
    └── (telemetry deposits; per-machine, not committed)
```

## Build

```bash
make all          # build modules + run all unit tests
make build        # compile only
make test         # run unit tests (requires built modules)
make cross-test   # run C2 cross-tests against Python fixtures (requires Python 3 + numpy)
make clean        # remove build artifacts
```

Requires **gfortran 13+** (any version supporting `-std=f2018 -pedantic`)
and **Python 3** with **NumPy** installed. The Makefile is OS-aware
(Windows cmd.exe + POSIX sh + macOS all supported via KB-038 pattern).

## Conventions

- **ASCII-only source** (KB-039): no Unicode in `.f90` files. Banner blocks use `=` × 57.
- **Lowercase loop indices** (KB-037): Fortran is case-insensitive, so `integer :: N, n`
  is a duplicate-declaration trap. Use `nlen` for sizes, lowercase for loop counters.
- **Equation-to-code mapping in module headers**: every source file opens
  with a doc-comment block showing how each canonical-tier equation maps
  to specific code lines.
- **No optimisation in v0.x reference builds**: `-O0 -ffp-contract=off`.

## C2 cross-test discipline

The cross-test driver compares Fortran outputs against `tests/fixtures/*.json`
(deterministic Python NumPy oracle) element-wise at
`tol = max(|ref|, floor) × ULP × sqrt(N)` per Higham 2002 §4.2 Option G.

Three cross-tests participate in the C2 gate:

| Test | Fixture | Cases | v0.1.0 worst residual |
|------|---------|-------|------------------------|
| gravity_gradient | gravity_gradient_inputs.json | 7 | 3.153e-16 (≈ 1 ULP) |
| multipole_potential | multipole_potential_inputs.json | 3 | 0.0 (exact) |
| symmetric_top | symmetric_top_trajectory.json | 201 sampled states | 3.110e-11 at step 12000 |

## License

Inherits from the parent project: Apache 2.0 for code, CC-BY-SA-4.0 for docs.
Per-routine source-header citations are part of the C3 acceptance criterion.
