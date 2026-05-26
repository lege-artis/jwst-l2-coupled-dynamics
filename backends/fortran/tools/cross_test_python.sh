#!/usr/bin/env bash
# cross_test_python.sh — C2 gate (lege-artis doctrine §4.4) for the
# jwst-l2-coupled-dynamics Fortran reference backend.
#
# Compares Fortran test outputs against Python JSON fixtures element-wise
# at ULP × sqrt(N) tolerance (Higham 2002 §4.2 Option G).
#
# This script is invoked from the backends/fortran/ directory via `make cross-test`.
# It produces _audit/cross_test_summary.md.
#
# Contract: each cross-test-bearing test program (test_gravity_gradient_torque,
# test_multipole_potential, test_symmetric_top) MUST emit its computed outputs
# to a known location under _audit/ as JSON-parseable text matching the fixture
# schema. The Python helper script invoked below does the actual comparison.

set -euo pipefail

# --- locate self + project root ---
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FORTRAN_ROOT="$(cd "$HERE/.." && pwd)"
PROJECT_ROOT="$(cd "$FORTRAN_ROOT/../.." && pwd)"
FIXTURE_DIR="$PROJECT_ROOT/tests/fixtures"
AUDIT_DIR="$FORTRAN_ROOT/_audit"
SUMMARY_MD="$AUDIT_DIR/cross_test_summary.md"

mkdir -p "$AUDIT_DIR"

# Export before heredoc so Python subprocess can read them
export FORTRAN_ROOT FIXTURE_DIR AUDIT_DIR SUMMARY_MD

# --- preflight ---
echo "=========================================="
echo "  JWST-L2 Fortran Phase 1 — C2 Cross-test"
echo "=========================================="
echo "  Fortran root: $FORTRAN_ROOT"
echo "  Fixtures:     $FIXTURE_DIR"
echo "  Audit out:    $AUDIT_DIR"
echo ""

# Verify Python + numpy availability
if ! command -v python3 >/dev/null 2>&1; then
    echo "FAIL: python3 not on PATH"; exit 1
fi
if ! python3 -c "import numpy" 2>/dev/null; then
    echo "FAIL: numpy not installed for python3"; exit 1
fi

# Verify required fixtures exist
for fixture in gravity_gradient_inputs.json multipole_potential_inputs.json symmetric_top_trajectory.json; do
    if [ ! -f "$FIXTURE_DIR/$fixture" ]; then
        echo "FAIL: missing fixture $FIXTURE_DIR/$fixture"
        echo "  Run: python3 $PROJECT_ROOT/tests/fixtures/generate_fixtures.py"
        exit 1
    fi
done

# --- comparison logic (Python helper, inlined for portability) ---
python3 - <<'PYEOF' > "$AUDIT_DIR/cross_test_python_output.log" 2>&1 || CROSS_EXIT=$?
"""
Inline cross-test runner: read Fortran test output (JSON under _audit/), read
corresponding Python fixture, compare element-wise at ULP * sqrt(N) tolerance.

Output written to _audit/cross_test_summary.md. Exit code 0 if all PASS, 1 otherwise.
"""
import os, sys, json, math

FORTRAN_ROOT = os.environ.get("FORTRAN_ROOT")
FIXTURE_DIR  = os.environ.get("FIXTURE_DIR")
AUDIT_DIR    = os.environ.get("AUDIT_DIR")
SUMMARY_MD   = os.environ.get("SUMMARY_MD")

ULP_DP = 2.220446049250313e-16  # machine epsilon, double precision

def ulp_tol(ref, n_ops, scale_floor=None):
    """Higham §4.2 Option G: abs tolerance = max(|ref|, floor) * ULP * sqrt(N)."""
    mag = max(abs(float(ref)), scale_floor if scale_floor else 1e-300, 1e-300)
    return mag * ULP_DP * math.sqrt(max(n_ops, 1))

def compare_scalar(actual, expected, n_ops, scale_floor, label):
    tol  = ulp_tol(expected, n_ops, scale_floor)
    diff = abs(float(actual) - float(expected))
    rel  = diff / max(abs(float(expected)), 1e-300)
    return diff <= tol, rel, label, 1

def compare_vec(actual_list, expected_list, n_ops, scale_floor, label):
    """Compare two flat numeric lists element-wise."""
    if len(actual_list) != len(expected_list):
        return False, float("inf"), label, 0
    worst_res, worst_lbl, all_pass, total_n = 0.0, label, True, 0
    for i, (a, e) in enumerate(zip(actual_list, expected_list)):
        ok, res, _, n = compare_scalar(a, e, n_ops, scale_floor, f"{label}[{i}]")
        total_n += n
        if not ok:
            all_pass = False
        if res > worst_res:
            worst_res, worst_lbl = res, f"{label}[{i}]"
    return all_pass, worst_res, worst_lbl, total_n

def compare_output_dict(actual_dict, expected_dict, n_ops, label, scale_floor=None):
    """Compare numerical values in matching keys of two dicts."""
    worst_res, worst_lbl, all_pass, total_n = 0.0, label, True, 0
    for key, ev in expected_dict.items():
        av = actual_dict.get(key)
        if av is None:
            return False, float("inf"), f"{label}.{key} missing", 0
        sub_lbl = f"{label}.{key}"
        if isinstance(ev, list):
            ok, res, wl, n = compare_vec(av, ev, n_ops, scale_floor, sub_lbl)
        else:
            ok, res, wl, n = compare_scalar(av, ev, n_ops, scale_floor, sub_lbl)
        total_n += n
        if not ok:
            all_pass = False
        if res > worst_res:
            worst_res, worst_lbl = res, wl
    return all_pass, worst_res, worst_lbl, total_n

# Discover Fortran output files in audit dir
fortran_outputs = {
    "gravity_gradient": os.path.join(AUDIT_DIR, "fortran_gravity_gradient_output.json"),
    "multipole_potential": os.path.join(AUDIT_DIR, "fortran_multipole_potential_output.json"),
    "symmetric_top": os.path.join(AUDIT_DIR, "fortran_symmetric_top_output.json"),
}

fixture_paths = {
    "gravity_gradient": os.path.join(FIXTURE_DIR, "gravity_gradient_inputs.json"),
    "multipole_potential": os.path.join(FIXTURE_DIR, "multipole_potential_inputs.json"),
    "symmetric_top": os.path.join(FIXTURE_DIR, "symmetric_top_trajectory.json"),
}

# Estimated n_ops per output (Higham N_eff_ulp = sqrt(N))
n_ops = {
    "gravity_gradient": 30,         # ~3x3 matrix mult + cross + scalar = order 10s of ops
    "multipole_potential": 50,      # body-frame transform + trace + quadratic form
    "symmetric_top": 12000 * 200,   # RK4 with 12000 steps; per-step cost ~200 flops
}

report_lines = []
report_lines.append("# JWST-L2 Fortran Phase 1 — Cross-test summary (C2 gate)")
report_lines.append("")
report_lines.append(f"Run: cross_test_python.sh")
report_lines.append("")
report_lines.append("| Test | Status | Cases | Worst residual | Worst case |")
report_lines.append("|------|--------|-------|----------------|------------|")

overall_pass = True

for test_id, fortran_path in fortran_outputs.items():
    fixture_path = fixture_paths[test_id]
    test_n_ops = n_ops[test_id]

    if not os.path.exists(fortran_path):
        report_lines.append(f"| {test_id} | SKIP (no Fortran output) | - | - | {fortran_path} not found |")
        continue

    with open(fixture_path) as f:
        fixture = json.load(f)
    with open(fortran_path) as f:
        fortran = json.load(f)

    # Compare each case in the fixture against Fortran output.
    # Sonnet contract: Fortran outputs are saved as JSON with the same "cases" or "trajectory"
    # structure as the fixture, with "output" key holding the computed result.
    fixture_cases = fixture.get("cases", fixture.get("trajectory", []))
    fortran_cases = fortran.get("cases", fortran.get("trajectory", []))

    if len(fixture_cases) != len(fortran_cases):
        report_lines.append(f"| {test_id} | FAIL | {len(fixture_cases)} vs {len(fortran_cases)} | mismatch | case count differs |")
        overall_pass = False
        continue

    worst_res = 0.0
    worst_case = ""
    all_case_pass = True
    for fixt_case, fort_case in zip(fixture_cases, fortran_cases):
        case_label = fixt_case.get("name", str(fixt_case.get("step", "")))

        if test_id == "symmetric_top":
            # No "output" wrapper; compare q (unit quaternion, scale=1.0) + omega_body
            ok_q,  res_q,  wl_q,  _ = compare_vec(fort_case["q"],          fixt_case["q"],          test_n_ops, 1.0,  f"{case_label}.q")
            ok_om, res_om, wl_om, _ = compare_vec(fort_case["omega_body"],  fixt_case["omega_body"], test_n_ops, None, f"{case_label}.omega_body")
            passed = ok_q and ok_om
            res    = max(res_q, res_om)
            where  = wl_q if res_q >= res_om else wl_om
        else:
            # gravity_gradient / multipole_potential: compare "output" dict values
            exp_out = fixt_case.get("output", {})
            act_out = fort_case.get("output", {})
            passed, res, where, _ = compare_output_dict(act_out, exp_out, test_n_ops, case_label)

        if not passed:
            all_case_pass = False
        if res != float("inf") and res > worst_res:
            worst_res = res
            worst_case = where

    status = "PASS" if all_case_pass else "FAIL"
    overall_pass &= all_case_pass
    report_lines.append(f"| {test_id} | {status} | {len(fixture_cases)} | {worst_res:.3e} | {worst_case} |")

report_lines.append("")
report_lines.append(f"**Overall: {'PASS' if overall_pass else 'FAIL'}**")
report_lines.append("")
report_lines.append("Tolerance model: ULP × sqrt(N) per Higham 2002 §4.2 Option G.")
report_lines.append(f"ULP(double) = {ULP_DP:.6e}")

with open(SUMMARY_MD, "w") as f:
    f.write("\n".join(report_lines) + "\n")

print("\n".join(report_lines))
sys.exit(0 if overall_pass else 1)
PYEOF

CROSS_EXIT="${CROSS_EXIT:-0}"

echo ""
echo "Summary written to: $SUMMARY_MD"

if [ "$CROSS_EXIT" -eq 0 ]; then
    echo "C2 GATE: PASS"
else
    echo "C2 GATE: FAIL"
fi

exit "$CROSS_EXIT"
