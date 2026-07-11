# Release notes — v0.3.1 (doc-patch: Ch 06 §6.5.3 erratum + elliptic-oracle regression guard + CS canonical tier complete + README refresh)

**Date:** 2026-07 (publish pending Pete gate action)
**Type:** documentation-patch release on top of v0.3.0. No dynamics, integrator, or Lyapunov code changes.

## Summary

v0.3.1 corrects an error in the published canonical-tier documentation, adds the
regression test whose absence allowed the error to ship, completes the public
Czech canonical tier, and refreshes the README that had lagged two releases.

## 1. Ch 06 §6.5.3 erratum (EN + CS)

The v0.1.0–v0.3.0 published boxed expressions for the Jacobi-elliptic
parametrisation of the torque-free asymmetric top were erroneous:

| Quantity | Published (wrong) | Corrected (verified) |
|---|---|---|
| Reduction prefactor | $(I_2-I_1)(I_3-I_2)/(I_1 I_2 I_3)$ | $(I_2-I_1)(I_3-I_2)/(I_1 I_2^2 I_3)$ |
| $L_y$ amplitude | $\alpha \approx 0.116$ (test config) | $\beta = \sqrt{I_2(2TI_3-\|\mathbf{L}\|^2)/(I_3-I_2)} \approx 2.007$ (~16× larger) |
| $\Omega^*$ | mixed-invariant form | $\sqrt{(I_3-I_2)(\|\mathbf{L}\|^2-2TI_1)/(I_1I_2I_3)}$ |
| $k^2$ | reciprocal-branch form | $(I_2-I_1)(2TI_3-\|\mathbf{L}\|^2)/[(I_3-I_2)(\|\mathbf{L}\|^2-2TI_1)]$ |
| Branch statement | $\|\mathbf{L}\|^2 < 2TI_2$ w/ dn on $L_z$ (self-contradictory) | $\|\mathbf{L}\|^2 > 2TI_2$ (polhode encircles $\hat{\mathbf{e}}_3$; matches the test configurations) |

Implied flip-period error of the published forms: 7–63% depending on
configuration. The corrected forms follow Landau-Lifshitz Vol I §37 and match a
`solve_ivp` (rtol 1e-12) integration of the exact Euler equations to ≲1e-12
relative in flip period and 8e-14 element-wise in $(L_x,L_y,L_z)(t)$ across the
§6.5.5 test configurations. A dated correction note is embedded at §6.5.3 in
both language mirrors.

**How it happened (honest record):** the flip-period fixture metadata described
an "elliptic $T_{LL}=4K(k)/\Omega^*$ cross-check" that was never implemented in
code; the numerical golden was self-consistent, so nothing exercised the boxed
formulas. Fourth boxed-formula drift incident (KB-043 class), first to reach a
public mirror. The §6.5.5 tolerance table now points to the real guard.

## 2. New: `tests/test_elliptic_oracle.py`

Six tests validating the documented closed form against tight-tolerance
numerics: flip period (≤1e-9), $L_y$ turning-point amplitude (≤1e-6),
element-wise trajectory match (≤1e-6; achieves ~1e-13), exact conservation
identities of the parametrisation (≤1e-13), and a two-point perturbation sweep
(≤1e-8). Doc boxes and test formulas must now move in lockstep — a drift fails
CI before it ships.

## 3. Czech canonical tier complete

`docs/canonical/cs/` now carries Ch 00–07 (this release adds 04, 05, 06 —
Eulerovy-Lagrangeovy rovnice, tělesové Eulerovy rovnice, analytické limity a
validační brány — the last including the §6.5.3 correction from day one).
Engineer-tier Ch 07 how-to available in EN + CS.

## 4. README refresh

The public README previously described the v0.1.0 state (listing the symplectic
integrator and Lyapunov spectrum as future work while both were in the published
tree). Status / Plan / documentation-tier tables now reflect v0.3.1 reality.

## Breaking changes

None. No numerical results change; only documentation and added tests.

## Known limitations / forward pointers

- The Master Test Set gate suite (19 oracle fixtures + reduced-ODE runner)
  remains in the source monorepo pending its convergence onto the internal
  test-orchestration framework (gate-wiring + naming corrections in progress);
  it will ship in a later release with honest scope labelling.
- L2 restricted three-body (CR3BP) extension: TG-8 oracles specified; targeted
  v0.4.x.

## Citation

Yamyang, P. (2026). lege-artis/jwst-l2-coupled-dynamics v0.3.1.
https://github.com/lege-artis/jwst-l2-coupled-dynamics (tag v0.3.1)
