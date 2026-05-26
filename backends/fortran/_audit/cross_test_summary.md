# JWST-L2 Fortran Phase 1 — Cross-test summary (C2 gate)

Run: cross_test_python.sh

| Test | Status | Cases | Worst residual | Worst case |
|------|--------|-------|----------------|------------|
| gravity_gradient | PASS | 7 | 3.153e-16 | case_04_asymmetric_body.tau_body[0] |
| multipole_potential | PASS | 3 | 0.000e+00 |  |
| symmetric_top | PASS | 201 | 3.110e-11 | 12000.q[1] |

**Overall: PASS**

Tolerance model: ULP × sqrt(N) per Higham 2002 §4.2 Option G.
ULP(double) = 2.220446e-16
