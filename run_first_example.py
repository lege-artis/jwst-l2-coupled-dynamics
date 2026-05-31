"""
run_first_example.py — runs the JWST-L2 first-cut scenario and writes pure
data outputs (no visualisation). NDJSON file per snapshot covers all the
quantities the strategic capture doc identified as worth showing.

Scenario (illustrative, first-cut):
  Body A: JWST-like composite, initially at the origin.
    Initial quaternion: 0.3 rad rotation about the inertial (1,1,0)/sqrt(2)
    axis (tilts body away from inertial alignment so the inertial-frame
    relative-position vector is NOT along a body principal axis — this
    activates the gravity-gradient torque).
    Initial angular velocity in body frame: omega = (0.02, 0.0, 0.08) rad/s.
    The non-zero perpendicular component (omega_x) combined with the parallel
    component (omega_z) drives free precession at the Euler frequency
    lambda = (I_zz - I_xx)/I_xx * omega_z ~ -0.027 rad/s ~ period 230 s.
  Body B: probe-like composite, initially at +x = 50 m from A, stationary
    translationally, with omega = (0.0, 0.03, 0.08) rad/s. Body B is much
    more elongated (inertia ratio ~5.8) so its precession period is shorter.
    Initial quaternion: identity.

Time horizon: 600 s with dt = 0.05 s -> 12000 steps. Snapshot every 20 steps
so the output file is ~600 records.

Output (NDJSON: one JSON object per line):
  {
    "t": float,                    -- time (s)
    "A": {
      "x": [3 floats],             -- position (inertial frame, m)
      "v": [3 floats],             -- velocity (inertial frame, m/s)
      "q": [4 floats],             -- quaternion (qw, qx, qy, qz)
      "R": [9 floats, row-major],  -- rotation matrix body -> inertial
      "omega_body": [3 floats],    -- angular velocity (body frame, rad/s)
      "omega_inertial": [3 floats],-- angular velocity (inertial frame, rad/s)
      "L_inertial": [3 floats],    -- spin angular momentum (inertial, kg m^2/s)
      "I_inertial": [9 floats],    -- inertia tensor in inertial frame
      "principal_axes_inertial": [9 floats],  -- 3 principal axes in inertial frame
      "principal_moments": [3 floats]         -- ascending eigenvalues of I
    },
    "B": {... same fields ...},
    "diagnostics": {
      "separation": float,
      "E_trans": float, "E_rot": float, "E_pot": float, "E_total": float,
      "L_total": [3 floats],
      "L_total_magnitude": float,
      "tau_gg_A_magnitude": float,  -- gravity-gradient torque on A (body frame, N m)
      "tau_gg_B_magnitude": float
    }
  }
"""
from __future__ import annotations
import json
import sys
import time
from pathlib import Path
import numpy as np

# Local imports (this file is meant to be invoked from the experiment directory)
sys.path.insert(0, str(Path(__file__).parent))
from geometries import make_jwst_like, make_probe_like
from dynamics import (
    BodyState, pack_state, unpack_state, rk4_step,
    q_to_matrix, inertia_inertial_frame, principal_axes_inertial,
    gravity_gradient_torque_body_frame,
    total_kinetic_energy, total_potential_energy, total_angular_momentum,
    back_reaction_force_one_body,
)


def main():
    # --- build bodies + inertia tensors in body frame ---
    body_a_composite = make_jwst_like()
    body_b_composite = make_probe_like()
    m_a = body_a_composite.total_mass
    m_b = body_b_composite.total_mass
    I_a_body = body_a_composite.inertia_at_com()
    I_b_body = body_b_composite.inertia_at_com()

    # --- initial conditions ---
    # Body A initial attitude: 0.3 rad rotation about inertial (1,1,0)/sqrt(2).
    # Quaternion: q = (cos(theta/2), sin(theta/2) * axis_unit)
    theta_a = 0.3
    axis_a = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
    q_a_init = np.concatenate([[np.cos(theta_a / 2.0)],
                                np.sin(theta_a / 2.0) * axis_a])
    initial_a = BodyState(
        x=[0.0, 0.0, 0.0],
        v=[0.0, 0.0, 0.0],
        q=q_a_init,
        omega=[0.02, 0.0, 0.08],          # perpendicular + parallel -> free precession
    )
    initial_b = BodyState(
        x=[50.0, 0.0, 0.0],               # 50 m along inertial +x
        v=[0.0, 0.0, 0.0],
        q=[1.0, 0.0, 0.0, 0.0],
        omega=[0.0, 0.03, 0.08],          # perpendicular (y) + parallel (z)
    )
    state = pack_state(initial_a, initial_b)

    # --- integration loop ---
    t = 0.0
    dt = 0.05         # 0.05 s timestep
    n_steps = 12000   # 600 s total horizon
    snapshot_every = 20

    outdir = Path(__file__).parent / "outputs"
    outdir.mkdir(exist_ok=True)
    outfile = outdir / "first_example_trajectory.ndjson"

    init_sep = float(np.linalg.norm(initial_b.x - initial_a.x))
    print(f"--- JWST-L2 first-cut scenario ---")
    print(f"  Body A (JWST-like):  m = {m_a:.2f} kg")
    print(f"  Body B (probe-like): m = {m_b:.2f} kg")
    print(f"  Initial separation:  {init_sep:.1f} m")
    print(f"  Initial omega_A (body): {initial_a.omega.tolist()} rad/s")
    print(f"  Initial omega_B (body): {initial_b.omega.tolist()} rad/s")
    print(f"  dt = {dt} s, n_steps = {n_steps} ({n_steps * dt:.1f} s horizon)")
    print(f"  Writing to: {outfile}")
    print()

    # --- baseline energy + angular momentum for conservation tracking ---
    Etk0, Erk0, _ = total_kinetic_energy(state, m_a, I_a_body, m_b, I_b_body)
    Ep0 = total_potential_energy(state, m_a, m_b, I_a_body, I_b_body)
    E0 = Etk0 + Erk0 + Ep0
    L0 = total_angular_momentum(state, m_a, I_a_body, m_b, I_b_body)
    print(f"  baseline E_trans = {Etk0:.6e} J")
    print(f"  baseline E_rot   = {Erk0:.6e} J")
    print(f"  baseline E_pot   = {Ep0:.6e} J")
    print(f"  baseline E_total = {E0:.6e} J")
    print(f"  baseline |L|     = {np.linalg.norm(L0):.6e} kg m^2/s")
    print()

    t_run0 = time.time()

    with open(outfile, "w") as f:
        for step in range(n_steps + 1):
            do_snapshot = (step % snapshot_every == 0)
            if do_snapshot:
                ba, bb = unpack_state(state)
                R_a = q_to_matrix(ba.q)
                R_b = q_to_matrix(bb.q)
                I_a_inertial = inertia_inertial_frame(I_a_body, ba.q)
                I_b_inertial = inertia_inertial_frame(I_b_body, bb.q)
                pmA_vals, pmA_vecs = principal_axes_inertial(I_a_body, ba.q)
                pmB_vals, pmB_vecs = principal_axes_inertial(I_b_body, bb.q)

                tau_a = gravity_gradient_torque_body_frame(
                    I_a_body, bb.x - ba.x, R_a, m_b)
                tau_b = gravity_gradient_torque_body_frame(
                    I_b_body, ba.x - bb.x, R_b, m_a)

                Etk, Erk, _ = total_kinetic_energy(state, m_a, I_a_body, m_b, I_b_body)
                Ep = total_potential_energy(state, m_a, m_b, I_a_body, I_b_body)
                L = total_angular_momentum(state, m_a, I_a_body, m_b, I_b_body)

                snap = {
                    "t": float(t),
                    "A": {
                        "x": ba.x.tolist(),
                        "v": ba.v.tolist(),
                        "q": ba.q.tolist(),
                        "R": R_a.flatten().tolist(),
                        "omega_body": ba.omega.tolist(),
                        "omega_inertial": (R_a @ ba.omega).tolist(),
                        "L_inertial": (R_a @ (I_a_body @ ba.omega)).tolist(),
                        "I_inertial": I_a_inertial.flatten().tolist(),
                        "principal_axes_inertial": pmA_vecs.flatten().tolist(),
                        "principal_moments": pmA_vals.tolist(),
                    },
                    "B": {
                        "x": bb.x.tolist(),
                        "v": bb.v.tolist(),
                        "q": bb.q.tolist(),
                        "R": R_b.flatten().tolist(),
                        "omega_body": bb.omega.tolist(),
                        "omega_inertial": (R_b @ bb.omega).tolist(),
                        "L_inertial": (R_b @ (I_b_body @ bb.omega)).tolist(),
                        "I_inertial": I_b_inertial.flatten().tolist(),
                        "principal_axes_inertial": pmB_vecs.flatten().tolist(),
                        "principal_moments": pmB_vals.tolist(),
                    },
                    "diagnostics": {
                        "separation": float(np.linalg.norm(bb.x - ba.x)),
                        "E_trans": float(Etk),
                        "E_rot": float(Erk),
                        "E_pot": float(Ep),
                        "E_total": float(Etk + Erk + Ep),
                        "L_total": L.tolist(),
                        "L_total_magnitude": float(np.linalg.norm(L)),
                        "tau_gg_A_magnitude": float(np.linalg.norm(tau_a)),
                        "tau_gg_B_magnitude": float(np.linalg.norm(tau_b)),
                    },
                }
                f.write(json.dumps(snap) + "\n")

            if step < n_steps:
                state = rk4_step(state, dt, m_a, I_a_body, m_b, I_b_body)
                t += dt

    t_run = time.time() - t_run0
    print(f"--- integration complete in {t_run:.2f} s ---")
    print(f"  snapshots written: {n_steps // snapshot_every + 1}")

    # --- final-state summary + conservation residuals ---
    ba, bb = unpack_state(state)
    Etk, Erk, _ = total_kinetic_energy(state, m_a, I_a_body, m_b, I_b_body)
    Ep = total_potential_energy(state, m_a, m_b, I_a_body, I_b_body)
    Ef = Etk + Erk + Ep
    Lf = total_angular_momentum(state, m_a, I_a_body, m_b, I_b_body)

    dE_rel = abs(Ef - E0) / abs(E0) if abs(E0) > 0 else abs(Ef - E0)
    dL_rel = np.linalg.norm(Lf - L0) / np.linalg.norm(L0) if np.linalg.norm(L0) > 0 else np.linalg.norm(Lf - L0)

    print()
    print(f"--- final state ---")
    print(f"  separation:      {np.linalg.norm(bb.x - ba.x):.3f} m  (initial {init_sep:.1f} m)")
    print(f"  omega_A (body):  {ba.omega}")
    print(f"  omega_B (body):  {bb.omega}")
    print()
    print(f"--- conservation residuals ---")
    print(f"  |dE / E0|: {dE_rel:.3e}  (target: ~RK4 truncation, ~ (dt)^4 = {dt**4:.3e})")
    print(f"  |dL / L0|: {dL_rel:.3e}")
    print()
    print(f"  Note: at L2 in the real problem, total angular momentum is NOT")
    print(f"  conserved because Sun + Earth exert external torques (in inertial")
    print(f"  frame). This first-cut prototype isolates only the A-B coupling")
    print(f"  (no Sun/Earth) so |dL/L0| above is a pure numerical-method check.")
    print()

    # --- significant rotational modes detection (quick frequency analysis) ---
    print(f"--- significant rotational modes (FFT-based, per-axis omega) ---")
    print(f"  This is a quick spectral peek to confirm precession/nutation")
    print(f"  modes are activated. Full lege-artis/fourier integration is")
    print(f"  the proper analysis pathway.")
    print()
    # read back the trajectory + FFT one axis
    omegas = []
    with open(outfile) as f:
        for line in f:
            snap = json.loads(line)
            omegas.append(snap["A"]["omega_body"])
    omegas = np.array(omegas)  # shape (n_snapshots, 3)
    # snapshot dt = dt * snapshot_every
    fs = 1.0 / (dt * snapshot_every)  # snapshot sampling frequency, Hz
    for axis_idx, axis_name in enumerate(("x", "y", "z")):
        signal = omegas[:, axis_idx] - omegas[:, axis_idx].mean()
        Xk = np.fft.rfft(signal)
        freqs = np.fft.rfftfreq(len(signal), d=1.0/fs)
        # exclude DC; find largest peak
        if len(freqs) > 1:
            pk_idx = 1 + int(np.argmax(np.abs(Xk[1:])))
            pk_freq = freqs[pk_idx]
            pk_mag = np.abs(Xk[pk_idx])
            print(f"  omega_A_body[{axis_name}]: peak freq = {pk_freq:.4f} Hz "
                  f"(period {1.0/pk_freq:.2f} s if non-zero), peak mag = {pk_mag:.3e}")

    print()
    print(f"--- done. Trajectory at: {outfile} ---")


if __name__ == "__main__":
    main()
