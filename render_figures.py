"""
render_figures.py — minimal sanity rendering of the JWST-L2 sweep outputs.

Decoupled from computation. Reads NDJSON trajectory files + SWEEP-SUMMARY.json
produced by testcases.py and emits a small set of figures showing the
data shape and key diagnostics. Per session reframing (2026-05-22), figures
are supplementary, not the headline deliverable — the structured-numerical
output is the value. This module is re-runnable any time the rendering needs
polish without re-running the dynamics.

Figures produced:
  Fig 1 — body B principal-moment anisotropy bar chart (per case)
  Fig 2 — omega-alignment_A_B time series (the cross-body coupling signature)
  Fig 3 — off-diagonal inertia magnitude evolution for body A (the tumble metric)
  Fig 4 — DFT spectrum of omega_A_body[x] for one representative case (C2)

All static PNG, deterministic, dpi 120, ≤ 300 KB each. No animation.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")               # headless backend
import matplotlib.pyplot as plt

OUTDIR = Path(__file__).parent / "outputs"
FIGDIR = OUTDIR / "figures"
FIGDIR.mkdir(exist_ok=True)


def load_summary():
    with open(OUTDIR / "SWEEP-SUMMARY.json") as f:
        return json.load(f)


def load_case(name: str) -> list[dict]:
    path = OUTDIR / f"case-{name}.ndjson"
    snaps = []
    with open(path) as f:
        for line in f:
            snaps.append(json.loads(line))
    return snaps


def fig1_anisotropy_bar(summary: dict) -> Path:
    cases = summary["cases"]
    names = [c["name"] for c in cases]
    aniso_B = [c["anisotropy_B"] for c in cases]
    aniso_A = [c["anisotropy_A"] for c in cases]

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(names))
    width = 0.4
    ax.bar(x - width/2, aniso_A, width, label="Body A (JWST-like)",
           color="#3a7ca5", edgecolor="black", linewidth=0.5)
    ax.bar(x + width/2, aniso_B, width, label="Body B (parametric cone)",
           color="#d9923f", edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Testcase")
    ax.set_ylabel("Inertia anisotropy I_max / I_min")
    ax.set_title("Per-case inertia anisotropy "
                 "(Body A constant; Body B varies with probe length)")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=30, ha="right", fontsize=9)
    ax.set_yscale("log")
    ax.grid(True, axis="y", alpha=0.3, which="both")
    ax.legend(loc="upper left")
    fig.tight_layout()
    outpath = FIGDIR / "fig1-anisotropy-per-case.png"
    fig.savefig(outpath, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return outpath


def fig2_omega_alignment(summary: dict) -> Path:
    fig, ax = plt.subplots(figsize=(10, 5))
    colours = plt.cm.tab10(np.linspace(0, 1, len(summary["cases"])))
    for case_summary, colour in zip(summary["cases"], colours):
        snaps = load_case(case_summary["name"])
        t = np.array([s["t"] for s in snaps])
        alignment = np.array([s["diagnostics"]["omega_alignment_A_B"]
                               for s in snaps])
        ax.plot(t, alignment, label=case_summary["name"], color=colour, lw=1.2)
    ax.axhline(0, color="black", lw=0.5, alpha=0.3)
    ax.set_xlabel("Simulation time (s)")
    ax.set_ylabel("cos(angle) between omega_A and omega_B (inertial frame)")
    ax.set_title("Cross-body angular-velocity alignment "
                 "(the coupling signature)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8, ncol=2)
    ax.set_ylim(-1.05, 1.05)
    fig.tight_layout()
    outpath = FIGDIR / "fig2-omega-alignment-timeseries.png"
    fig.savefig(outpath, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return outpath


def fig3_off_diagonal_body_A(summary: dict) -> Path:
    fig, ax = plt.subplots(figsize=(10, 5))
    colours = plt.cm.tab10(np.linspace(0, 1, len(summary["cases"])))
    for case_summary, colour in zip(summary["cases"], colours):
        snaps = load_case(case_summary["name"])
        t = np.array([s["t"] for s in snaps])
        off_A = np.array([s["A"]["off_diag_magnitude"] for s in snaps])
        ax.plot(t, off_A, label=case_summary["name"], color=colour, lw=1.2)
    ax.set_xlabel("Simulation time (s)")
    ax.set_ylabel("RMS off-diagonal magnitude of I_A (inertial, kg.m^2)")
    ax.set_title("Body A — off-diagonal inertia-tensor magnitude in inertial frame "
                 "(tumbling signature)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=8, ncol=2)
    fig.tight_layout()
    outpath = FIGDIR / "fig3-off-diagonal-body-A.png"
    fig.savefig(outpath, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return outpath


def fig4_dft_omega_A_x(case_name: str = "C2-medium-x") -> Path:
    snaps = load_case(case_name)
    omega_x = np.array([s["A"]["omega_body"][0] for s in snaps])
    t = np.array([s["t"] for s in snaps])
    dt_snap = t[1] - t[0]
    fs = 1.0 / dt_snap
    signal = omega_x - omega_x.mean()
    Xk = np.fft.rfft(signal)
    freqs = np.fft.rfftfreq(len(signal), d=dt_snap)
    mag = np.abs(Xk)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7))
    # time-domain (zoom to first 600 s for readability)
    mask_t = t <= 600
    ax1.plot(t[mask_t], omega_x[mask_t], color="#3a7ca5", lw=1.2)
    ax1.set_xlabel("Simulation time (s)")
    ax1.set_ylabel("omega_A_body[x] (rad/s)")
    ax1.set_title(f"Time series of omega_A_body[x] — case {case_name}")
    ax1.grid(True, alpha=0.3)
    # frequency-domain
    ax2.semilogy(freqs[1:], mag[1:], color="#d9923f", lw=1.0)
    pk_idx = 1 + int(np.argmax(mag[1:]))
    ax2.axvline(freqs[pk_idx], color="red", lw=0.7, linestyle="--",
                label=f"observed peak {freqs[pk_idx]:.4f} Hz "
                      f"(period {1.0/freqs[pk_idx]:.1f} s)")
    # theoretical Euler frequency line
    # body A: pmA = [15384, 23322, 23322]; ω_z = 0.08
    lambda_theory = (15384.0 - 23322.0) / 23322.0 * 0.08
    f_theory = abs(lambda_theory) / (2.0 * np.pi)
    ax2.axvline(f_theory, color="green", lw=0.7, linestyle=":",
                label=f"theoretical Euler frequency {f_theory:.4f} Hz "
                      f"(period {1.0/f_theory:.1f} s)")
    ax2.set_xlabel("Frequency (Hz)")
    ax2.set_ylabel("|DFT(omega_A_body[x])| (rad/s)")
    ax2.set_title(f"Power spectrum of omega_A_body[x] — case {case_name}")
    ax2.grid(True, alpha=0.3, which="both")
    ax2.legend(loc="best", fontsize=9)
    ax2.set_xlim(0, min(0.05, freqs[-1]))    # focus on low-frequency band
    fig.tight_layout()
    outpath = FIGDIR / "fig4-dft-spectrum-case-C2.png"
    fig.savefig(outpath, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return outpath


def main():
    summary = load_summary()
    print(f"--- rendering minimal sanity figures from {len(summary['cases'])} testcases ---")
    paths = [
        fig1_anisotropy_bar(summary),
        fig2_omega_alignment(summary),
        fig3_off_diagonal_body_A(summary),
        fig4_dft_omega_A_x(),
    ]
    print(f"--- rendered {len(paths)} figures to {FIGDIR} ---")
    for p in paths:
        size_kb = p.stat().st_size / 1024
        print(f"  {p.name}: {size_kb:.1f} KB")


if __name__ == "__main__":
    main()
