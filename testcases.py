"""
testcases.py — parametric sweep of the JWST-L2 first-cut dynamics.

Six testcases demonstrating non-trivial coupling between two "almost free"
bodies, designed to reveal that the gravity-gradient torque + inertia-tensor
geometry produce measurable coupling effects even at small force magnitudes.

Sweep structure:
  - Cases C1-C3: probe length varies (1.0 m, 5.0 m, 15.0 m), approach
    direction fixed along inertial +x. Body B inertia anisotropy sweeps
    from 3.25 (short rod) through 78 to 703 (long pencil).
  - Cases C4-C6: probe length fixed at 5.0 m, approach direction varies
    (+y axis, +z axis, diagonal (1,1,1)/sqrt(3)). Same anisotropy; different
    relative geometry between A's principal axes and the inter-body vector.

The two sweep axes isolate two distinct coupling mechanisms:
  - Length sweep -> how anisotropy magnitude affects coupling response
  - Direction sweep -> how approach-vs-principal-axis alignment affects
                       gravity-gradient torque direction

All cases use:
  - Body A: make_jwst_like() with tilted initial orientation
  - Body B: parametric cone probe at the indicated length
  - Initial separation 30 m (closer than the first-example 50 m, to amplify
    gravity-gradient torque visibility)
  - Integration horizon 6000 s, dt = 0.05 s, snapshot every 100 steps
    (long enough for GG-induced drift to surface on top of free precession)

Per-case output: outputs/case-CN.ndjson + per-case summary in
outputs/SWEEP-SUMMARY.json.
"""
from __future__ import annotations
import json
import sys
import time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from geometries import make_jwst_like, make_parametric_cone_probe
from dynamics import (
    BodyState, pack_state, unpack_state, rk4_step,
    q_to_matrix, inertia_inertial_frame, principal_axes_inertial,
    gravity_gradient_torque_body_frame,
    total_kinetic_energy, total_potential_energy, total_angular_momentum,
)


# --- testcase parameter table ---

TESTCASES = [
    # name,           probe_L,  approach_direction_unit_vector
    ("C1-short-rod-x",   1.0,  np.array([1.0, 0.0, 0.0])),
    ("C2-medium-x",      5.0,  np.array([1.0, 0.0, 0.0])),
    ("C3-long-pencil-x", 15.0, np.array([1.0, 0.0, 0.0])),
    ("C4-medium-y",      5.0,  np.array([0.0, 1.0, 0.0])),
    ("C5-medium-z",      5.0,  np.array([0.0, 0.0, 1.0])),
    ("C6-medium-diag",   5.0,  np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)),
]

SEPARATION_M = 30.0           # initial inter-body separation
INTEGRATION_HORIZON_S = 1000  # 16.7 min — sufficient for free precession + GG-drift onset
DT_S = 0.1                    # RK4 timestep (period > 100 s -> dt = 1/1000 of period -> plenty)
SNAPSHOT_EVERY = 20           # one snapshot per 2 s of simulated time -> 501 snapshots/case
OUTDIR = Path(__file__).parent / "outputs"


# --- diagnostic helpers ---

def off_diagonal_magnitude(I_inertial_3x3: np.ndarray) -> float:
    """RMS magnitude of the off-diagonal elements of a 3x3 symmetric tensor.
       For a diagonal-in-body-frame inertia tensor, this quantity grows from
       zero when the body is body-axis-aligned to a maximum when fully tumbled."""
    off = np.array([
        I_inertial_3x3[0, 1], I_inertial_3x3[0, 2], I_inertial_3x3[1, 2]
    ])
    return float(np.sqrt(np.mean(off * off)))


def cross_body_omega_alignment(omega_A_inertial: np.ndarray,
                                omega_B_inertial: np.ndarray) -> float:
    """Cosine of the angle between the two bodies' angular velocity vectors
       in the inertial frame. Close to +/-1 = aligned (potentially resonant).
       Close to 0 = perpendicular (no rotational coupling expected)."""
    mag_a = np.linalg.norm(omega_A_inertial)
    mag_b = np.linalg.norm(omega_B_inertial)
    if mag_a < 1e-12 or mag_b < 1e-12:
        return 0.0
    return float(np.dot(omega_A_inertial, omega_B_inertial) / (mag_a * mag_b))


# --- one testcase runner ---

def run_one_testcase(name: str, probe_L: float, approach_unit: np.ndarray):
    print(f"\n=== testcase {name}: probe L = {probe_L} m, approach = {approach_unit.tolist()} ===")

    # bodies + inertia tensors
    body_a_composite = make_jwst_like()
    body_b_composite = make_parametric_cone_probe(length=probe_L)
    m_a = body_a_composite.total_mass
    m_b = body_b_composite.total_mass
    I_a_body = body_a_composite.inertia_at_com()
    I_b_body = body_b_composite.inertia_at_com()
    pmA = np.linalg.eigvalsh(I_a_body)
    pmB = np.linalg.eigvalsh(I_b_body)

    # initial conditions
    # Body A: tilted attitude as in run_first_example (activates GG torque)
    theta_a = 0.3
    axis_a = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0)
    q_a_init = np.concatenate([[np.cos(theta_a / 2.0)],
                                np.sin(theta_a / 2.0) * axis_a])
    initial_a = BodyState(
        x=[0.0, 0.0, 0.0],
        v=[0.0, 0.0, 0.0],
        q=q_a_init,
        omega=[0.02, 0.0, 0.08],   # mixed perp + parallel -> free precession
    )
    # Body B: positioned along the approach direction at SEPARATION_M
    initial_b = BodyState(
        x=(SEPARATION_M * approach_unit).tolist(),
        v=[0.0, 0.0, 0.0],
        q=[1.0, 0.0, 0.0, 0.0],
        omega=[0.0, 0.03, 0.08],
    )
    state = pack_state(initial_a, initial_b)

    # baselines
    Etk0, Erk0, _ = total_kinetic_energy(state, m_a, I_a_body, m_b, I_b_body)
    Ep0 = total_potential_energy(state, m_a, m_b)
    E0 = Etk0 + Erk0 + Ep0
    L0 = total_angular_momentum(state, m_a, I_a_body, m_b, I_b_body)

    # outfile
    outfile = OUTDIR / f"case-{name}.ndjson"

    # integration loop
    t = 0.0
    n_steps = int(INTEGRATION_HORIZON_S / DT_S)
    snapshot_count = 0
    max_off_diag_A = 0.0
    max_off_diag_B = 0.0
    max_tau_gg_A = 0.0
    max_tau_gg_B = 0.0
    initial_alignment = None
    final_alignment = None
    t_run0 = time.time()

    with open(outfile, "w") as f:
        for step in range(n_steps + 1):
            do_snapshot = (step % SNAPSHOT_EVERY == 0)
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
                Ep = total_potential_energy(state, m_a, m_b)
                L = total_angular_momentum(state, m_a, I_a_body, m_b, I_b_body)

                offA = off_diagonal_magnitude(I_a_inertial)
                offB = off_diagonal_magnitude(I_b_inertial)
                omega_A_inertial = R_a @ ba.omega
                omega_B_inertial = R_b @ bb.omega
                alignment = cross_body_omega_alignment(omega_A_inertial, omega_B_inertial)
                if initial_alignment is None:
                    initial_alignment = alignment
                final_alignment = alignment

                tau_a_mag = float(np.linalg.norm(tau_a))
                tau_b_mag = float(np.linalg.norm(tau_b))
                max_off_diag_A = max(max_off_diag_A, offA)
                max_off_diag_B = max(max_off_diag_B, offB)
                max_tau_gg_A = max(max_tau_gg_A, tau_a_mag)
                max_tau_gg_B = max(max_tau_gg_B, tau_b_mag)

                snap = {
                    "t": float(t),
                    "A": {
                        "x": ba.x.tolist(),
                        "v": ba.v.tolist(),
                        "q": ba.q.tolist(),
                        "omega_body": ba.omega.tolist(),
                        "omega_inertial": omega_A_inertial.tolist(),
                        "I_inertial": I_a_inertial.flatten().tolist(),
                        "off_diag_magnitude": offA,
                        "principal_moments": pmA_vals.tolist(),
                    },
                    "B": {
                        "x": bb.x.tolist(),
                        "v": bb.v.tolist(),
                        "q": bb.q.tolist(),
                        "omega_body": bb.omega.tolist(),
                        "omega_inertial": omega_B_inertial.tolist(),
                        "I_inertial": I_b_inertial.flatten().tolist(),
                        "off_diag_magnitude": offB,
                        "principal_moments": pmB_vals.tolist(),
                    },
                    "diagnostics": {
                        "separation": float(np.linalg.norm(bb.x - ba.x)),
                        "E_trans": float(Etk),
                        "E_rot": float(Erk),
                        "E_pot": float(Ep),
                        "E_total": float(Etk + Erk + Ep),
                        "L_total_magnitude": float(np.linalg.norm(L)),
                        "tau_gg_A_magnitude": tau_a_mag,
                        "tau_gg_B_magnitude": tau_b_mag,
                        "omega_alignment_A_B": alignment,
                    },
                }
                f.write(json.dumps(snap) + "\n")
                snapshot_count += 1

            if step < n_steps:
                state = rk4_step(state, DT_S, m_a, I_a_body, m_b, I_b_body)
                t += DT_S

    t_run = time.time() - t_run0

    # final state diagnostic
    ba, bb = unpack_state(state)
    Etk_f, Erk_f, _ = total_kinetic_energy(state, m_a, I_a_body, m_b, I_b_body)
    Ep_f = total_potential_energy(state, m_a, m_b)
    Ef = Etk_f + Erk_f + Ep_f
    Lf = total_angular_momentum(state, m_a, I_a_body, m_b, I_b_body)
    dE_rel = abs(Ef - E0) / abs(E0) if abs(E0) > 0 else abs(Ef - E0)
    dL_rel = (np.linalg.norm(Lf - L0) / np.linalg.norm(L0)
              if np.linalg.norm(L0) > 0 else np.linalg.norm(Lf - L0))

    # quick FFT of omega_A_body[x] to extract the precession peak
    omegas = []
    with open(outfile) as fh:
        for line in fh:
            omegas.append(json.loads(line)["A"]["omega_body"])
    omegas = np.array(omegas)
    fs = 1.0 / (DT_S * SNAPSHOT_EVERY)
    signal_x = omegas[:, 0] - omegas[:, 0].mean()
    Xk = np.fft.rfft(signal_x)
    freqs = np.fft.rfftfreq(len(signal_x), d=1.0/fs)
    pk_idx = 1 + int(np.argmax(np.abs(Xk[1:]))) if len(freqs) > 1 else 0
    pk_freq = float(freqs[pk_idx])
    pk_period = float(1.0 / pk_freq) if pk_freq > 0 else 0.0

    # theoretical Euler frequency for body A (axisymmetric body, principal moments pmA)
    omega_z_init = 0.08
    if pmA[0] != pmA[-1]:
        # axial moment is the smallest one (oblate A) -> I_zz - I_xx negative
        I_axial = pmA[0]
        I_perp = pmA[-1]
        lambda_theory = (I_axial - I_perp) / I_perp * omega_z_init
        period_theory = float(2.0 * np.pi / abs(lambda_theory)) if lambda_theory != 0 else 0.0
    else:
        period_theory = 0.0

    summary = {
        "name": name,
        "probe_L": probe_L,
        "approach_direction": approach_unit.tolist(),
        "separation_initial_m": SEPARATION_M,
        "integration_horizon_s": INTEGRATION_HORIZON_S,
        "snapshots": snapshot_count,
        "wall_time_s": round(t_run, 2),
        "m_A_kg": m_a,
        "m_B_kg": m_b,
        "principal_moments_A": pmA.tolist(),
        "principal_moments_B": pmB.tolist(),
        "anisotropy_A": float(pmA[-1] / pmA[0]),
        "anisotropy_B": float(pmB[-1] / pmB[0]),
        "fft_peak_period_omega_A_x": pk_period,
        "fft_peak_period_theory_axisymmetric_A": period_theory,
        "max_off_diag_A": max_off_diag_A,
        "max_off_diag_B": max_off_diag_B,
        "max_tau_gg_A_Nm": max_tau_gg_A,
        "max_tau_gg_B_Nm": max_tau_gg_B,
        "initial_omega_alignment": initial_alignment,
        "final_omega_alignment": final_alignment,
        "dE_over_E0": dE_rel,
        "dL_over_L0": dL_rel,
    }
    print(f"  wall time: {t_run:.2f} s; snapshots: {snapshot_count}")
    print(f"  |dE/E0|: {dE_rel:.3e}; |dL/L0|: {dL_rel:.3e}")
    print(f"  body B anisotropy: {summary['anisotropy_B']:.2f}")
    print(f"  FFT peak period for omega_A_x: {pk_period:.2f} s (theory: {period_theory:.2f} s)")
    print(f"  max gravity-gradient torque on A: {max_tau_gg_A:.3e} N.m")
    print(f"  max off-diagonal magnitude in I_A_inertial: {max_off_diag_A:.3e} kg.m^2")
    print(f"  omega-alignment A-B: initial {initial_alignment:.3f} -> final {final_alignment:.3f}")
    return summary


# --- driver ---

def main():
    OUTDIR.mkdir(exist_ok=True)
    print(f"--- parametric sweep over {len(TESTCASES)} testcases ---")
    print(f"  Each case: {INTEGRATION_HORIZON_S}-s horizon at dt={DT_S}, "
          f"separation {SEPARATION_M} m")
    print(f"  Output: {OUTDIR}/case-CN.ndjson + SWEEP-SUMMARY.json")

    sweep_t0 = time.time()
    summaries = []
    for (name, probe_L, approach) in TESTCASES:
        summaries.append(run_one_testcase(name, probe_L, approach))
    sweep_t = time.time() - sweep_t0

    summary_file = OUTDIR / "SWEEP-SUMMARY.json"
    with open(summary_file, "w") as f:
        json.dump({
            "sweep_total_wall_time_s": round(sweep_t, 2),
            "n_testcases": len(summaries),
            "cases": summaries,
        }, f, indent=2)
    print(f"\n--- sweep complete in {sweep_t:.2f} s; summary at {summary_file} ---")


if __name__ == "__main__":
    main()
