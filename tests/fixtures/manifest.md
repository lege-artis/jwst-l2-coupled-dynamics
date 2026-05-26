# JWST-L2 Fortran Phase 1 — Cross-test fixture manifest

Generated: 2026-05-24T09:55:07.680013
Python: 3.10.12
NumPy: 2.2.6
G_NEWTON: 6.6743e-11 (m^3 / kg / s^2; CODATA 2018)

## Files

### gravity_gradient_inputs.json
- spec_id: gravity_gradient_inputs.v0.1
- 7 test cases for D4 cross-test
- Each case has (I_body, r_rel_inertial, R_body_to_inertial, m_other) input + tau_body output
- Cases cover: libration equilibrium, off-equilibrium tilt, probe attitude, asymmetric body, near-spherical edge case

### symmetric_top_trajectory.json
- spec_id: symmetric_top_trajectory.v0.1
- Torque-free Euler precession of axisymmetric JWST-like body for D8 cross-test
- Setup: I = diag(2540, 2540, 3870), omega_initial = (1e-3, 0, 1e-2) rad/s, dt = 0.05 s, n_steps = 12000 (window 600 s)
- Body B set to negligible mass + far away to suppress gravity-gradient torque
- Trajectory sampled every 60 steps -> ~200 state vectors of body A
- Predicted Euler precession frequency lambda = (I_par-I_perp)/I_perp * omega_z = 5.236220e-03 rad/s
- Predicted libration period T = 1199.9467 s
- Achieved |dE/E_0| = 8.3956e-42
- Achieved |dL/L_0| = 1.5572e-15

### multipole_potential_inputs.json
- spec_id: multipole_potential_inputs.v0.1
- 3 test cases for D3 cross-test
- Each case has (m_A, m_B, rho, I_A, q_A, I_B, q_B) input + (V_mono, V_tidal_A, V_tidal_B, V_total) output
- Cases cover: far-field, near-field, point-mass Kepler limit

## How the Fortran cross-test consumes these

Run the Fortran D3/D4/D8 test programs against the same inputs. Compare element-
wise outputs to the JSON values at ULP * sqrt(N) tolerance per doctrine §4.4 C2.
The cross-test runner script `backends/fortran/tools/cross_test_python.sh` is
the canonical harness; it loads each fixture, runs the corresponding Fortran
test, and emits PASS/FAIL telemetry to `backends/fortran/_audit/cross_test_summary.md`.

## Regeneration

This fixture set is reproducible from a single command:
    python tests/fixtures/generate_fixtures.py

The output is deterministic (no random seeds; all inputs are explicit). If a
fixture needs to be regenerated (e.g., after a dynamics.py revision), commit
the new JSON files in the same commit that touches dynamics.py.

## Authority

- Canonical source for gravity_gradient_torque: dynamics.py::gravity_gradient_torque_body_frame
- Canonical source for V_mono + V_quad: Chapter 03 §3.2 + §3.4.3 (this script
  implements the boxed formulas directly; dynamics.py::total_potential_energy
  only includes V_mono, so V_quad is generated from the canonical-tier formula
  in v_quadrupole_tidal_one_body() above)
- Canonical source for RK4 trajectory: dynamics.py::rk4_step + state_derivative
