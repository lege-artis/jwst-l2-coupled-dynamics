"""
generate_fixtures.py — Phase 1 cross-test oracle generation.

Produces JSON fixtures consumed by the Fortran lege-artis Phase 1 cross-tests
(D4 gravity-gradient torque + D8 symmetric-top trajectory + D3 multipole potential).
Per SONNET-BRIEF.md §4: this script is the Pete-side / Opus-side authority for the
C2 oracle data.

Run from the project root (experiments/jwst-l2-first-cut/):
    python tests/fixtures/generate_fixtures.py

Outputs:
    tests/fixtures/gravity_gradient_inputs.json
    tests/fixtures/symmetric_top_trajectory.json
    tests/fixtures/multipole_potential_inputs.json
    tests/fixtures/manifest.md

All numerics are IEEE-754 double precision. JSON uses repr-equivalent float
serialisation so that round-trip reading recovers the bit-identical values
(important: the C2 cross-test tolerance is ULP × sqrt(N), so any fixture-side
precision loss eats into the gate budget).
"""
from __future__ import annotations
import json
import os
import sys
import datetime
import numpy as np

# Import the project's first-cut modules
HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, PROJECT_ROOT)

import dynamics
import geometries

OUT_DIR = HERE


# --- helpers ---

def to_jsonable(x):
    """Convert numpy arrays to nested lists for JSON serialisation, preserving
       full double-precision via repr."""
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (np.floating, float)):
        return float(x)
    if isinstance(x, (np.integer, int)):
        return int(x)
    return x


def q_from_axis_angle(axis, angle):
    """Unit quaternion (qw, qx, qy, qz) for rotation by `angle` about `axis`."""
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    c = np.cos(angle / 2.0)
    s = np.sin(angle / 2.0)
    return np.array([c, s * axis[0], s * axis[1], s * axis[2]])


# --- multipole potential closed form (the canonical formula from Ch 03 §3.4.3) ---

def v_monopole(m_a, m_b, rho):
    """Monopole-monopole: V_mono = -G m_A m_B / |rho|."""
    return -dynamics.G_NEWTON * m_a * m_b / np.linalg.norm(rho)


def v_quadrupole_tidal_one_body(m_other, rho, I_body, q):
    """Quadrupole tidal contribution for ONE body (per Ch 03 §3.4.3 boxed):
       V_tidal,X = -(G m_other / (2 |rho|^3)) * [tr(I_X) - 3 rho_hat_body . I_X . rho_hat_body]"""
    r_mag = np.linalg.norm(rho)
    R = dynamics.q_to_matrix(q)
    rho_hat = rho / r_mag
    rho_hat_body = R.T @ rho_hat
    tr_I = np.trace(I_body)
    quadratic = rho_hat_body @ I_body @ rho_hat_body
    return -(dynamics.G_NEWTON * m_other) / (2.0 * r_mag ** 3) * (tr_I - 3.0 * quadratic)


def v_total_truncated(m_a, m_b, rho, I_a, q_a, I_b, q_b):
    """Total quadrupole-truncated mutual potential."""
    v0 = v_monopole(m_a, m_b, rho)
    vA = v_quadrupole_tidal_one_body(m_b, rho, I_a, q_a)
    # For body B: rho points from A to B; from B's perspective it's -rho
    vB = v_quadrupole_tidal_one_body(m_a, -rho, I_b, q_b)
    return {"v_mono": v0, "v_tidal_A": vA, "v_tidal_B": vB,
            "v_total": v0 + vA + vB}


# --- fixture 1: gravity-gradient torque inputs/outputs ---

def make_gravity_gradient_fixture():
    """5 representative torque-computation cases covering equilibrium, off-equilibrium,
       axisymmetric vs asymmetric bodies, near-zero corner cases."""
    cases = []

    # Case 1: JWST-like body, libration equilibrium (rho along symmetry axis)
    jwst = geometries.make_jwst_like()
    I_a = jwst.inertia_at_com()
    q_a = np.array([1.0, 0.0, 0.0, 0.0])  # identity
    rho = np.array([0.0, 0.0, 1.5e9])     # along inertial z, matching body z
    R = dynamics.q_to_matrix(q_a)
    tau = dynamics.gravity_gradient_torque_body_frame(I_a, rho, R, 5.972e24)
    cases.append({
        "name": "case_01_jwst_libration_eq",
        "description": "JWST-like body at hypothetical L2-like configuration, rho_hat along body z-axis (libration equilibrium). Expected tau = 0.",
        "inputs": {
            "I_body": to_jsonable(I_a),
            "r_rel_inertial": to_jsonable(rho),
            "R_body_to_inertial": to_jsonable(R),
            "m_other": 5.972e24,
        },
        "output": {"tau_body": to_jsonable(tau)},
    })

    # Case 2: JWST-like body, off-equilibrium (10 deg tilt about y-axis)
    q_a = q_from_axis_angle([0.0, 1.0, 0.0], np.radians(10.0))
    R = dynamics.q_to_matrix(q_a)
    tau = dynamics.gravity_gradient_torque_body_frame(I_a, rho, R, 5.972e24)
    cases.append({
        "name": "case_02_jwst_off_eq_10deg_y",
        "description": "JWST-like body tilted 10 deg about inertial y; rho still along inertial z; oblate -> destabilising torque about body y expected.",
        "inputs": {
            "I_body": to_jsonable(I_a),
            "r_rel_inertial": to_jsonable(rho),
            "R_body_to_inertial": to_jsonable(R),
            "m_other": 5.972e24,
        },
        "output": {"tau_body": to_jsonable(tau)},
    })

    # Case 3: Probe-like body at generic attitude
    probe = geometries.make_probe_like()
    I_b = probe.inertia_at_com()
    q_b = q_from_axis_angle([1.0, 1.0, 1.0], np.radians(30.0))
    rho_probe = np.array([1e6, 2e6, 5e8])  # generic inertial direction
    R_b = dynamics.q_to_matrix(q_b)
    tau_probe = dynamics.gravity_gradient_torque_body_frame(I_b, rho_probe, R_b, 5.972e24)
    cases.append({
        "name": "case_03_probe_generic_attitude",
        "description": "Probe-like body, q ~ 30 deg about (1,1,1) axis; rho in arbitrary inertial direction.",
        "inputs": {
            "I_body": to_jsonable(I_b),
            "r_rel_inertial": to_jsonable(rho_probe),
            "R_body_to_inertial": to_jsonable(R_b),
            "m_other": 5.972e24,
        },
        "output": {"tau_body": to_jsonable(tau_probe)},
    })

    # Case 4: Fully asymmetric body, generic configuration
    I_asym = np.diag([1000.0, 2000.0, 5000.0])
    q_asym = q_from_axis_angle([0.3, 0.7, 0.5], np.radians(45.0))
    rho_asym = np.array([3e8, -2e8, 4e8])
    R_asym = dynamics.q_to_matrix(q_asym)
    tau_asym = dynamics.gravity_gradient_torque_body_frame(I_asym, rho_asym, R_asym, 5.972e24)
    cases.append({
        "name": "case_04_asymmetric_body",
        "description": "Fully asymmetric inertia diag(1000, 2000, 5000); generic rho and attitude. Exercises all three principal-axis distinctness.",
        "inputs": {
            "I_body": to_jsonable(I_asym),
            "r_rel_inertial": to_jsonable(rho_asym),
            "R_body_to_inertial": to_jsonable(R_asym),
            "m_other": 5.972e24,
        },
        "output": {"tau_body": to_jsonable(tau_asym)},
    })

    # Case 5: Near-spherical body, expect near-zero torque
    I_near_sph = np.diag([2540.0, 2541.0, 2540.5])  # nearly equal
    rho_near = np.array([0.0, 0.0, 1.5e9])
    q_near = np.array([1.0, 0.0, 0.0, 0.0])
    R_near = dynamics.q_to_matrix(q_near)
    tau_near = dynamics.gravity_gradient_torque_body_frame(I_near_sph, rho_near, R_near, 5.972e24)
    cases.append({
        "name": "case_05_near_spherical",
        "description": "Near-spherical body (anisotropy ~ 5e-4). Expected tau magnitude many orders below the principal-moment-difference scale.",
        "inputs": {
            "I_body": to_jsonable(I_near_sph),
            "r_rel_inertial": to_jsonable(rho_near),
            "R_body_to_inertial": to_jsonable(R_near),
            "m_other": 5.972e24,
        },
        "output": {"tau_body": to_jsonable(tau_near)},
    })

    # Cases 6 + 7: Oblate reference body (per OQ-FORT-6 Path C resolution 2026-05-24).
    # Pairs with the make_jwst_like prolate composite to demonstrate both regimes of
    # the gravity-gradient stability table from Ch 03 Sec.3.7.2.
    oblate = geometries.make_oblate_reference_body()
    I_oblate = oblate.inertia_at_com()  # diag(2540, 2540, 3870) by construction

    # Case 6: oblate body at libration equilibrium (rho along body z-axis).
    # Expected torque is exactly zero (the equilibrium configuration).
    q_oblate_eq = np.array([1.0, 0.0, 0.0, 0.0])
    rho_oblate = np.array([0.0, 0.0, 1.5e9])
    R_oblate_eq = dynamics.q_to_matrix(q_oblate_eq)
    tau_oblate_eq = dynamics.gravity_gradient_torque_body_frame(I_oblate, rho_oblate, R_oblate_eq, 5.972e24)
    cases.append({
        "name": "case_06_oblate_libration_eq",
        "description": "Oblate reference body diag(2540,2540,3870) at libration equilibrium (rho_hat along body z). Expected tau = 0. This is the OBLATE counterpart of case_01_jwst_libration_eq; the JWST-like composite is prolate by parallel-axis spread, while this synthetic body is oblate by construction.",
        "inputs": {
            "I_body": to_jsonable(I_oblate),
            "r_rel_inertial": to_jsonable(rho_oblate),
            "R_body_to_inertial": to_jsonable(R_oblate_eq),
            "m_other": 5.972e24,
        },
        "output": {"tau_body": to_jsonable(tau_oblate_eq)},
    })

    # Case 7: oblate body 10 deg off-equilibrium about y. With I_par > I_perp (oblate),
    # the gravity-gradient torque should DESTABILISE (push the symmetry axis away from
    # radial), so tau_y should be nonzero with the destabilising sign.
    q_oblate_off = q_from_axis_angle([0.0, 1.0, 0.0], np.radians(10.0))
    R_oblate_off = dynamics.q_to_matrix(q_oblate_off)
    tau_oblate_off = dynamics.gravity_gradient_torque_body_frame(I_oblate, rho_oblate, R_oblate_off, 5.972e24)
    cases.append({
        "name": "case_07_oblate_off_eq_10deg_y",
        "description": "Oblate reference body tilted 10 deg about inertial y; rho still along inertial z. Oblate body -> destabilising torque about body y expected, sign opposite to case_02_jwst_off_eq_10deg_y (the prolate JWST-like would be restoring at this attitude).",
        "inputs": {
            "I_body": to_jsonable(I_oblate),
            "r_rel_inertial": to_jsonable(rho_oblate),
            "R_body_to_inertial": to_jsonable(R_oblate_off),
            "m_other": 5.972e24,
        },
        "output": {"tau_body": to_jsonable(tau_oblate_off)},
    })

    return {
        "spec_id": "gravity_gradient_inputs.v0.1",
        "generated_at": datetime.datetime.now().isoformat(),
        "python_version": sys.version,
        "numpy_version": np.__version__,
        "G_NEWTON": dynamics.G_NEWTON,
        "n_cases": len(cases),
        "cases": cases,
    }


# --- fixture 2: torque-free symmetric-top trajectory (D8 C2 gate) ---

def make_symmetric_top_trajectory():
    """Torque-free Euler precession of an axisymmetric body. Set body B "absent"
       (m_B negligible and far away) so the gravity-gradient torque vanishes.
       Integrate body A for 600 s at dt = 0.05 s and save the full trajectory."""
    # JWST-like inertia (axisymmetric)
    I_a = np.diag([2540.0, 2540.0, 3870.0])
    m_a = 2450.0

    # Initial body-frame angular velocity: omega_perp + omega_z
    omega_perp = 1.0e-3   # rad/s
    omega_z = 1.0e-2      # rad/s

    body_a = dynamics.BodyState(
        x=np.array([0.0, 0.0, 0.0]),
        v=np.array([0.0, 0.0, 0.0]),
        q=np.array([1.0, 0.0, 0.0, 0.0]),
        omega=np.array([omega_perp, 0.0, omega_z]),
    )

    # Body B "absent": very small mass + very far away.
    I_b = np.eye(3) * 1e-20
    m_b = 1e-20  # effectively zero
    body_b = dynamics.BodyState(
        x=np.array([1e15, 0.0, 0.0]),
        v=np.array([0.0, 0.0, 0.0]),
        q=np.array([1.0, 0.0, 0.0, 0.0]),
        omega=np.array([0.0, 0.0, 0.0]),
    )

    dt = 0.05  # s
    n_steps = 12000  # 600 s window

    state = dynamics.pack_state(body_a, body_b)

    # Record initial diagnostics
    E_trans0, E_rot0, E_total0 = dynamics.total_kinetic_energy(state, m_a, I_a, m_b, I_b)
    V0 = dynamics.total_potential_energy(state, m_a, m_b)
    L0 = dynamics.total_angular_momentum(state, m_a, I_a, m_b, I_b)
    L0_mag = float(np.linalg.norm(L0))

    # Save downsampled trajectory: every 60th step (10 Hz sampling -> 200 samples)
    # to keep JSON manageable; the full trajectory's 12001 states with 13 doubles
    # per state per body = ~2 MB JSON. Downsampling to every 60th gives ~33 kB.
    sample_every = 60
    trajectory = []
    for step in range(n_steps + 1):
        if step % sample_every == 0:
            body_a_now, body_b_now = dynamics.unpack_state(state)
            # Record body A state only (body B is the absent one)
            trajectory.append({
                "step": step,
                "t": step * dt,
                "x": to_jsonable(body_a_now.x),
                "v": to_jsonable(body_a_now.v),
                "q": to_jsonable(body_a_now.q),
                "omega_body": to_jsonable(body_a_now.omega),
            })
        if step < n_steps:
            state = dynamics.rk4_step(state, dt, m_a, I_a, m_b, I_b)

    # Final-state diagnostics
    body_a_final, body_b_final = dynamics.unpack_state(state)
    E_trans_f, E_rot_f, E_total_f = dynamics.total_kinetic_energy(state, m_a, I_a, m_b, I_b)
    V_f = dynamics.total_potential_energy(state, m_a, m_b)
    L_f = dynamics.total_angular_momentum(state, m_a, I_a, m_b, I_b)

    # Predicted Euler precession frequency (closed-form, LL §35)
    lambda_pred = (I_a[2, 2] - I_a[0, 0]) / I_a[0, 0] * omega_z

    return {
        "spec_id": "symmetric_top_trajectory.v0.1",
        "generated_at": datetime.datetime.now().isoformat(),
        "python_version": sys.version,
        "numpy_version": np.__version__,
        "setup": {
            "I_a_body": to_jsonable(I_a),
            "m_a": m_a,
            "omega_initial_body": [omega_perp, 0.0, omega_z],
            "q_initial": [1.0, 0.0, 0.0, 0.0],
            "I_b_body_negligible": to_jsonable(I_b),
            "m_b_negligible": m_b,
            "rho_initial": [1e15, 0.0, 0.0],
            "dt": dt,
            "n_steps": n_steps,
            "sample_every": sample_every,
        },
        "diagnostics": {
            "E_total_initial": E_total0 + V0,
            "E_total_final": E_total_f + V_f,
            "dE_over_E0": (E_total_f + V_f - E_total0 - V0) / (E_total0 + V0) if E_total0 + V0 != 0 else None,
            "L_initial_inertial": to_jsonable(L0),
            "L_final_inertial": to_jsonable(L_f),
            "dL_over_L0": float(np.linalg.norm(L_f - L0) / L0_mag),
            "lambda_predicted_rad_per_s": lambda_pred,
            "lambda_T_predicted_s": 2.0 * np.pi / abs(lambda_pred),
        },
        "trajectory": trajectory,
    }


# --- fixture 3: multipole potential evaluations ---

def make_multipole_potential_fixture():
    """Evaluate V_mono + V_quad_tidal_A + V_quad_tidal_B for several
       (rho, I_A, q_A, I_B, q_B) configurations. The Fortran D3 implementation
       must reproduce these element-wise."""
    cases = []

    jwst = geometries.make_jwst_like()
    probe = geometries.make_probe_like()
    I_a = jwst.inertia_at_com()
    I_b = probe.inertia_at_com()
    m_a = jwst.total_mass
    m_b = probe.total_mass

    # Case 1: bodies far apart (>> body sizes), generic orientations
    rho = np.array([0.0, 0.0, 1e6])
    q_a = np.array([1.0, 0.0, 0.0, 0.0])
    q_b = q_from_axis_angle([1.0, 0.0, 0.0], np.radians(45.0))
    v = v_total_truncated(m_a, m_b, rho, I_a, q_a, I_b, q_b)
    cases.append({
        "name": "case_01_far_apart_generic",
        "description": "JWST + probe at 1000 km separation; far field; quadrupole correction is tiny.",
        "inputs": {
            "m_a": m_a, "m_b": m_b,
            "rho": to_jsonable(rho),
            "I_a_body": to_jsonable(I_a),
            "q_a": to_jsonable(q_a),
            "I_b_body": to_jsonable(I_b),
            "q_b": to_jsonable(q_b),
        },
        "output": {k: float(val) for k, val in v.items()},
    })

    # Case 2: bodies closer but still well-separated
    rho = np.array([0.0, 0.0, 1000.0])
    v = v_total_truncated(m_a, m_b, rho, I_a, q_a, I_b, q_b)
    cases.append({
        "name": "case_02_close_1km",
        "description": "JWST + probe at 1 km separation; quadrupole correction is large fraction of monopole.",
        "inputs": {
            "m_a": m_a, "m_b": m_b,
            "rho": to_jsonable(rho),
            "I_a_body": to_jsonable(I_a),
            "q_a": to_jsonable(q_a),
            "I_b_body": to_jsonable(I_b),
            "q_b": to_jsonable(q_b),
        },
        "output": {k: float(val) for k, val in v.items()},
    })

    # Case 3: point masses (I_A = I_B = 0) -- pure Kepler limit
    rho = np.array([0.0, 0.0, 1.5e9])
    q_a = np.array([1.0, 0.0, 0.0, 0.0])
    q_b = np.array([1.0, 0.0, 0.0, 0.0])
    v = v_total_truncated(m_a, m_b, rho, np.zeros((3, 3)), q_a, np.zeros((3, 3)), q_b)
    cases.append({
        "name": "case_03_point_masses_kepler_limit",
        "description": "I_A = I_B = 0 (point masses); V_total should equal V_mono exactly, V_tidal = 0.",
        "inputs": {
            "m_a": m_a, "m_b": m_b,
            "rho": to_jsonable(rho),
            "I_a_body": to_jsonable(np.zeros((3, 3))),
            "q_a": to_jsonable(q_a),
            "I_b_body": to_jsonable(np.zeros((3, 3))),
            "q_b": to_jsonable(q_b),
        },
        "output": {k: float(val) for k, val in v.items()},
    })

    return {
        "spec_id": "multipole_potential_inputs.v0.1",
        "generated_at": datetime.datetime.now().isoformat(),
        "python_version": sys.version,
        "numpy_version": np.__version__,
        "G_NEWTON": dynamics.G_NEWTON,
        "n_cases": len(cases),
        "cases": cases,
    }


# --- manifest writer ---

def write_manifest(grav_data, traj_data, potential_data):
    manifest = f"""# JWST-L2 Fortran Phase 1 — Cross-test fixture manifest

Generated: {datetime.datetime.now().isoformat()}
Python: {sys.version.split()[0]}
NumPy: {np.__version__}
G_NEWTON: {dynamics.G_NEWTON} (m^3 / kg / s^2; CODATA 2018)

## Files

### gravity_gradient_inputs.json
- spec_id: {grav_data['spec_id']}
- {grav_data['n_cases']} test cases for D4 cross-test
- Each case has (I_body, r_rel_inertial, R_body_to_inertial, m_other) input + tau_body output
- Cases cover: libration equilibrium, off-equilibrium tilt, probe attitude, asymmetric body, near-spherical edge case

### symmetric_top_trajectory.json
- spec_id: {traj_data['spec_id']}
- Torque-free Euler precession of axisymmetric JWST-like body for D8 cross-test
- Setup: I = diag(2540, 2540, 3870), omega_initial = (1e-3, 0, 1e-2) rad/s, dt = 0.05 s, n_steps = 12000 (window 600 s)
- Body B set to negligible mass + far away to suppress gravity-gradient torque
- Trajectory sampled every 60 steps -> ~200 state vectors of body A
- Predicted Euler precession frequency lambda = (I_par-I_perp)/I_perp * omega_z = {traj_data['diagnostics']['lambda_predicted_rad_per_s']:.6e} rad/s
- Predicted libration period T = {traj_data['diagnostics']['lambda_T_predicted_s']:.4f} s
- Achieved |dE/E_0| = {traj_data['diagnostics']['dE_over_E0']:.4e}
- Achieved |dL/L_0| = {traj_data['diagnostics']['dL_over_L0']:.4e}

### multipole_potential_inputs.json
- spec_id: {potential_data['spec_id']}
- {potential_data['n_cases']} test cases for D3 cross-test
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
"""
    with open(os.path.join(OUT_DIR, "manifest.md"), "w") as f:
        f.write(manifest)


# --- main driver ---

def main():
    print(f"Generating fixtures in {OUT_DIR}/")

    print("  - gravity_gradient_inputs.json ...")
    grav = make_gravity_gradient_fixture()
    with open(os.path.join(OUT_DIR, "gravity_gradient_inputs.json"), "w") as f:
        json.dump(grav, f, indent=2)
    print(f"    {grav['n_cases']} cases written")

    print("  - symmetric_top_trajectory.json ...")
    traj = make_symmetric_top_trajectory()
    with open(os.path.join(OUT_DIR, "symmetric_top_trajectory.json"), "w") as f:
        json.dump(traj, f, indent=2)
    print(f"    {len(traj['trajectory'])} sampled states (out of {traj['setup']['n_steps']} integrated)")
    print(f"    |dE/E_0| = {traj['diagnostics']['dE_over_E0']:.4e}")
    print(f"    |dL/L_0| = {traj['diagnostics']['dL_over_L0']:.4e}")

    print("  - multipole_potential_inputs.json ...")
    potential = make_multipole_potential_fixture()
    with open(os.path.join(OUT_DIR, "multipole_potential_inputs.json"), "w") as f:
        json.dump(potential, f, indent=2)
    print(f"    {potential['n_cases']} cases written")

    print("  - manifest.md ...")
    write_manifest(grav, traj, potential)

    print("\nAll fixtures generated.")


if __name__ == "__main__":
    main()
