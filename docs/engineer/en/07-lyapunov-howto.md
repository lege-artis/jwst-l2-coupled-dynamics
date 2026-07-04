# Engineer-tier how-to — computing Lyapunov spectra

> **Audience.** A practitioner who wants to *run* the Lyapunov-spectrum tools in `lyapunov/` against a trajectory and read off the answer — not derive the algorithm. For the first-principles derivation (BGGS, the R^14 embedding, the integrability criterion), see the canonical tier: `docs/canonical/en/07-lyapunov-spectrum.md`.
>
> **Status.** v0.1 (covers Phase 2 single-vector + Phase 3 full-spectrum tooling). First chapter of the `docs/engineer/en/` tier.
>
> **TL;DR.** Want the sign of chaos (is λ_max > 0)? Use `lyapunov_spectrum(..., n_vectors=1)`. Want the *full* 12-dimensional spectrum + Hamiltonian pairing? Use `lyapunov_spectrum_deflated(...)` (cheap) or `lyapunov_spectrum_liealg(...)` (manifold-native). Don't use the flat R^14 driver at `n_vectors=14` for a full spectrum — it injects two spurious contracting modes (see §5).

---

## §1 When to reach for which tool

| You want… | Call | Cost |
|---|---|---|
| λ_max only (sign of chaos) | `lyapunov.spectrum.lyapunov_spectrum(..., n_vectors=1)` | 1 trajectory + 1 tangent vector |
| Full 12-dim spectrum, cheap | `lyapunov.cocycle_deflated.lyapunov_spectrum_deflated(...)` | 1 trajectory + 12 tangent vectors, R^14 FD + deflation |
| Full 12-dim spectrum, manifold-native | `lyapunov.cocycle_liealg.lyapunov_spectrum_liealg(...)` | same count; so(3)×R^3 cocycle |
| Calibrate the pipeline against known physics | `testbeds.lyapunov_three_scenario` (a/b/c) | three reference runs |
| Reproduce the Phase 3 full-spectrum gates | `testbeds.lyapunov_full_spectrum` | Tier A/B runners |

All four drivers share the same call shape and return the same `dict` (see §4). The integrator `step_fn` is normally `integrators.symplectic_se3_variational.se3_variational_step` — the symplectic KDK step that conserves angular momentum to ULP, so the reference trajectory is clean.

## §2 Quick start — λ_max for the Dzhanibekov scenario

```python
import numpy as np
from dynamics import BodyState, pack_state
from integrators.symplectic_se3_variational import se3_variational_step
from lyapunov.spectrum import lyapunov_spectrum
from lyapunov.reference import dzhanibekov_characteristic_frequency

# Asymmetric body near the intermediate-axis saddle, Earth-mass primary nearby.
I_a = np.diag([100.0, 200.0, 400.0])      # I1 < I2 < I3 (required for Dzhanibekov)
m_a = 500.0
omega0 = np.array([1e-3, 1e-2, 1e-3])     # near e2 (intermediate axis)
m_b = 5.97e24
I_b = np.eye(3) * (2/5 * m_b * 6.371e6**2)
rho = 3.0e6

body_a = BodyState(np.zeros(3), np.zeros(3), np.array([1.,0,0,0]), omega0)
body_b = BodyState(np.array([rho,0,0]), np.zeros(3), np.array([1.,0,0,0]), np.zeros(3))
state0 = pack_state(body_a, body_b)

result = lyapunov_spectrum(
    state0, se3_variational_step, dt=50.0, n_steps=100_000, n_vectors=1,
    m_a=m_a, I_a_body=I_a, m_b=m_b, I_b_body=I_b,
)

lam_max = result['exponents'][0]
sigma   = dzhanibekov_characteristic_frequency(100., 200., 400., omega0[1])
print(f"lambda_max = {lam_max:.3e} s^-1   (sigma = {sigma:.3e})")
# lambda_max ~ 1.5e-5 > 0  -> weakly chaotic (gravity-gradient broke integrability)
```

A positive `lam_max` means deterministic chaos; `lam_max ≈ 0` with bounded motion means regular/integrable. Here ~1.5e-5 > 0 confirms the gravity-gradient coupling makes the near-saddle motion weakly chaotic.

## §3 Full 12-dimensional spectrum

For the whole spectrum (and the Hamiltonian-pairing check), use a 12-dimensional cocycle. Two equivalent entry points — both take the **same args as `lyapunov_spectrum`** and default `n_vectors=12`:

```python
from lyapunov.cocycle_deflated import lyapunov_spectrum_deflated   # Tier A
from lyapunov.cocycle_liealg  import lyapunov_spectrum_liealg      # Tier B

res = lyapunov_spectrum_deflated(
    state0, se3_variational_step, dt=50.0, n_steps=4_000,
    m_a=m_a, I_a_body=I_a, m_b=m_b, I_b_body=I_b,
)
print(res['exponents'])        # 12 values, descending
print(res['sum_exponents'])    # should be << 0.14 (near 0 for a Hamiltonian flow)
```

- **Tier A** (`_deflated`) keeps the R^14 finite-difference cocycle but projects each tangent vector orthogonal to the two quaternion-radial directions (which co-evolve with the trajectory) before each Gram-Schmidt. Cheapest path to a clean 12-dim spectrum.
- **Tier B** (`_liealg`) propagates tangents in the body-frame Lie algebra so(3)×R^3, which is 12-dimensional *by construction* — the quaternion-norm direction never appears. Same evaluation count (1 + 2·n_vectors integrator calls per step). Use this when you want the manifold-native answer with no post-hoc projection.

Both give a spectrum that sums to ≈ 0 and supports Hamiltonian pairing (§6).

## §4 Reading the result dict

Every driver returns the same keys:

| Key | Meaning |
|---|---|
| `exponents` | `n_vectors` Lyapunov exponents, descending (s⁻¹) |
| `sum_exponents` | sum of computed exponents (full-spectrum drivers; ≈ 0 expected) |
| `tau_qr_used` | production QR reorthonormalisation cadence (s), autotuned |
| `lambda_max_pilot` | λ_max estimate from the cheap pilot run |
| `n_qr` | number of QR steps performed |
| `running_exponents` | `[(t, exponents), …]` snapshots — plot these to see convergence |
| `T_total` | total simulated time in the production run (s) |
| `final_state` | 26-component state at the end |

The QR cadence is autotuned: a short pilot run estimates λ_max, then the production cadence is set to ~0.1 e-foldings per epoch (clamped). You normally don't set `tau_qr` by hand — pass `None` (the default) and let it tune.

## §5 The one thing not to do

**Do not call the flat `lyapunov_spectrum(..., n_vectors=14)` and read it as the full spectrum.** The flat R^14 chart carries two non-dynamical quaternion-norm directions (one per body) that the integrator's per-step normalisation contracts hard. They show up as two large-negative exponents (≈ −0.064, −0.071 s⁻¹) with no positive partners, so the R^14 spectrum sums to ≈ −0.14 instead of 0 and **fails Hamiltonian pairing**. The single-vector λ_max is unaffected (it tracks the most-expanding *dynamical* direction), which is why `n_vectors=1` is fine — but for the full spectrum always use the 12-dim deflated or Lie-algebraic driver. Full diagnosis: canonical Ch 07 §7.7 + §7.10; empirical record `_audit/phase3-prep/FULL-SPECTRUM-FINDING-2026-06-07.md`.

## §6 Gates — checking your spectrum is sane

Helpers in `lyapunov.reference` and `lyapunov.pairing` turn a spectrum into pass/fail checks:

```python
from lyapunov.pairing import spectrum_sum_gate, spurious_modes_absent, full_spectrum_gates
from lyapunov.reference import hamiltonian_pairing_check

spectrum_sum_gate(res['exponents'])        # |sum| < 1e-2  (Hamiltonian volume preservation)
spurious_modes_absent(res['exponents'])    # no exponent < -1e-2 (the norm modes are gone)
hamiltonian_pairing_check(res['exponents'])# {'pairing_pass': True, 'n_near_zero': ...}
full_spectrum_gates(res['exponents'])      # all of the above bundled
```

For the single-vector λ_max gates against known physics, `lyapunov.reference` has the scenario gates:
`lambda_max_gate_scenario_a` (integrable: λ_max < 1e-3), `lambda_max_gate_scenario_b` (libration: |λ_max| < ω_lib), `lambda_max_gate_scenario_c` (Dzhanibekov: λ_max > 0).

## §7 Gotchas

- **Integrable cases converge slowly.** For a regular/integrable trajectory the finite-time λ_max decays to 0 only as O(log T / T), so you'll see e.g. ~2e-4 at T=20,000 s rather than ~0. This is **intrinsic to the quasi-periodic dynamics, not the embedding** — the Lie-algebraic driver shows the same slow decay (Phase 3 confirmed this; canonical Ch 07 §7.3.2/§7.10). Don't set an integrable-case threshold tighter than ~1e-3, and judge "zero" against the relevant dynamical frequency (σ, ω_lib), not against 0.
- **Plot `running_exponents`** before trusting a single final number — if it hasn't plateaued, integrate longer.
- **Wall time** scales with n_vectors. A 12-vector run over ~4,000 steps is tens of seconds; budget accordingly (use `production_run=False` in the testbed runners for a quick smoke check).
- **Scenario (c) chaos is weak** (λ_c ~ 1.5e-5; divergence time ~1/λ ≈ 6.6e4 s) — near-separatrix KAM-regime, not strong stochasticity. Small-but-positive is the expected, correct answer.

## §8 Where to go next

- **Why** the algorithm is what it is, and the integrability criterion: canonical `docs/canonical/en/07-lyapunov-spectrum.md` (CS: `docs/canonical/cs/07-ljapunovovo-spektrum.md`).
- **The three calibration scenarios** in full: `testbeds/lyapunov_three_scenario.py`.
- **Phase 3 full-spectrum runners + gates**: `testbeds/lyapunov_full_spectrum.py`, `tests/test_lyapunov_phase3.py`.

---

**End 07-lyapunov-howto.md (engineer-tier v0.1).**
