# JWST-L2 Phase 2 — integrator landscape (deep-research output)

**Status:** v0.1 — deep-research landscape report; informs the Phase 2 Sonnet brief at `_handoffs/phase-2-symplectic-lie-group/SONNET-BRIEF.md`.  
**Author:** Opus (deep-research workflow), 2026-06-01.  
**Trigger:** Project-owner directive after WORKSTREAM A finalize close — "next shot" → JWST-L2 Phase 2 integrator landscape.  
**License:** CC-BY-SA-4.0 (companion doc).  
**Sibling artefacts:**
- `_handoffs/phase-2-symplectic-lie-group/SONNET-BRIEF.md` v0.1 (locked 2026-05-26 PM) — the implementation contract this landscape informs
- `docs/canonical/en/04-euler-lagrange-and-euler-poincare.md` — Lagrangian formulation
- `docs/canonical/en/06-analytical-limits-and-validation-gates.md` §6.4.3 closed-form Lagrange-top + §6.5.3 Jacobi-elliptic Dzhanibekov — analytical oracles this landscape positions against
- `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md` — sibling testbed-and-edge-case scoping

---

## §0. How to read this report

This is the *research-context layer* for Phase 2. It triangulates the published prior art on symplectic Lie-group integrators for coupled rigid-body dynamics, on Lyapunov-spectrum algorithms, on the Dzhanibekov flip in dynamical-systems literature, and on the JWST-specific orbital + attitude context. The existing Sonnet brief at `_handoffs/phase-2-symplectic-lie-group/SONNET-BRIEF.md` is the *implementation contract*; this landscape provides the *why this algorithm, why this tolerance, why this testbed* layer underneath the contract.

Two reading orders work:

1. **For the Sonnet executor:** read this landscape first (§1 → §5), then the SONNET-BRIEF.md, then apply the brief.
2. **For project owner reviewing Phase 2:** skim §6 (open questions) and §7 (recommended brief amendments) first, decide on the project-owner-decision points, then dive into §2 if any integrator-choice question opens up.

§1 frames Phase 2 against Phase 1 inheritance. §2 surveys integrator-choice prior art with a comparison table. §3 covers libration testbed design. §4 covers Lyapunov-spectrum methodology. §5 covers C2 cross-test tolerance calibration specific to symplectic integrators. §6 surfaces open questions. §7 recommends brief amendments.

---

## §1. Framing + Phase 1 inheritance

### §1.1 What Phase 1 established

Phase 1 (closed `jwst-l2-fortran-v0.1.0-phase-1-green` 2026-05-25) shipped:
- **Python reference layer** — 12-DOF coupled translational + rotational dynamics; quaternion-based attitude; classical Runge-Kutta 4 integrator; gravity-gradient torque per chapter 03 mutual gravitational potential; 43-test pytest suite; conservation diagnostics `|dE/E₀| ≈ 6e-12` and `|dL/L₀| ≈ 2e-10` at integration windows of order 1000 seconds with timestep 0.01 s.
- **Fortran reference backend** — 5 modules (`jwst_l2_{constants,quaternion,geometries,potential,torque,jacobian,diagnostics,dynamics}`) + 13 unit tests + 3 C2 cross-tests (D4 gravity-gradient torque, D5 multipole potential, D8 symmetric-top trajectory). C2 tolerance is **Higham §4.2 Option G ULP×√N** (Higham, *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM 2002 §4.2 Eq. 4.6 — the ULP-bounded probabilistic error model for floating-point sums). Worst-case relative error at green tag: ~3.1e-11 on the D8 symmetric-top trajectory at step 12000.
- **Canonical-tier docs Ch 00–06** — covering configuration manifold on SE(3) × SE(3), kinetic energy in intrinsic form, mutual gravitational potential (Cartesian + spherical-harmonic, cross-validated at quadrupole order), Euler-Lagrange + Euler-Poincaré reduction (Approach A in Euler-angle chart + Approach B via Lie-group reduction, cross-validated), body-frame Euler equations, and analytical-limits + validation gates. Ch 06 specifically includes §6.4.3 closed-form Lagrange-top + gravity-gradient libration solution via Podolský 2018 §8.3 effective-potential method, and §6.5.3 closed-form Jacobi-elliptic-function parametrisation of torque-free asymmetric top per Landau-Lifshitz Vol. I §37.

### §1.2 What Phase 2 needs to add

The Phase 1 RK4 integrator does its job at the trajectory scale (≤1000 s) — energy drift is at the IEEE-754 ULP-budget ceiling, angular-momentum drift only 2 orders above. The pedagogical and technical case for Phase 2 is the **long-time regime** (≥10⁴ s for libration-period observations, ≥10⁶ s for Lyapunov-spectrum convergence near Dzhanibekov-flip parameter regions). RK4 drifts in both energy and momentum at rate ∝ Δt^4 per step accumulated, so over 10⁸ time steps you accumulate enough drift to obscure the dynamics the long-time regime is supposed to reveal. A symplectic Lie-group integrator removes that long-time-drift floor.

Phase 2 thus adds three capabilities, per the Sonnet brief §1 table:

| Capability | What | Why |
|---|---|---|
| **Symplectic Lie-group variational integrator** | Discrete-Euler-Poincaré integrator on SO(3) × SO(3) (with SE(3) translations split off into the gravitational two-body) | Provable structure preservation; bounded energy oscillation under backward-error-analysis arguments instead of secular drift |
| **Libration testbed** | Long-time integration of a near-equilibrium configuration validated against the Ch 06 §6.4.3 closed-form Lagrange-top + gravity-gradient libration period | Confirms the integrator's qualitative + quantitative agreement with a *closed-form analytical oracle* at integration lengths where RK4 starts to drift visibly |
| **Lyapunov-spectrum methods** | Benettin-Galgani-Giorgilli-Strelcyn 1980 reorthogonalisation method (QR cadence) applied to SO(3) × SO(3) tangent-space cocycle | Quantitative chaos-characterisation near the Dzhanibekov intermediate-axis-instability regime |

---

## §2. Integrator choices

The literature offers three broad families of structure-preserving integrators for rigid-body dynamics on Lie groups. The names below are the canonical labels; each cites the foundational reference.

### §2.1 Variational integrators on Lie groups (Marsden-West / Lee-Leok-McClamroch family)

**Foundational references:**
- **Marsden, Pekarsky, Shkoller (1999)**, "Discrete Euler-Poincaré and Lie-Poisson equations," *Nonlinearity* **12**, 1647–1662. The original construction: a G-invariant discrete Lagrangian on a Lie group G, with discrete Hamilton's principle producing the discrete Euler-Poincaré equations, and the corresponding Lie-Poisson reduction.
- **Marsden, Pekarsky, Shkoller (2000)**, "Symmetry reduction of discrete Lagrangian mechanics on Lie groups," *J. Geom. Phys.* **36**, 139–150. Companion paper covering the symmetry-reduction theory.
- **Marsden, West (2001)**, "Discrete mechanics and variational integrators," *Acta Numerica* **10**, 357–514. **The canonical survey of variational integrators** — referenced extensively across the field, including the JWST-relevant Lee-Leok-McClamroch line of work below.

**Application to rigid-body attitude dynamics:**
- **Lee, Leok, McClamroch (2005)**, "A Lie group variational integrator for the attitude dynamics of a rigid body with applications to the 3D pendulum," *Proc. IEEE Conf. Control Applications*, pp. 962–967. The construction for a single rigid body on SO(3).
- **Lee, Leok, McClamroch (2007)**, "Lie group variational integrators for the full body problem in orbital mechanics," *Celestial Mechanics and Dynamical Astronomy* **98**, 121–144. DOI 10.1007/s10569-007-9073-x. **This paper is the direct prior-art reference for the JWST-L2 full body problem.** The abstract states: *"Equations of motion, referred to as full body models, are developed to describe the dynamics of rigid bodies acting under their mutual gravitational potential ... discrete equations of motion, referred to as a Lie group variational integrator, provide a geometrically exact and numerically efficient computational method for simulating full body dynamics in orbital mechanics; they are symplectic and momentum preserving, and they exhibit good energy behavior for exponentially long time periods ... efficient in only requiring a single evaluation of the gravity forces and moments per time step. The Lie group variational integrator also preserves the group structure without the use of local charts, reprojection, or constraints. Computational results are given for the dynamics of two rigid dumbbell bodies acting under their mutual gravity ..."* As-cited 98 times per SpringerLink as of fetch date.

The Lee-Leok-McClamroch (2007) construction is the **single most relevant prior art** for Phase 2: it applies exactly the variational-integrator-on-SE(3)×SE(3) machinery to the same physics (two rigid bodies under mutual gravitational potential) Phase 1's Python ref-impl already simulates. The Sonnet brief should anchor here.

**Other variational-tier references worth mentioning in `refs.bib`:**
- **Iserles, Munthe-Kaas, Nørsett, Zanna (2000)**, "Lie-group methods," *Acta Numerica* **9**, 215–365. The authoritative survey of Lie-group integration methods, including RKMK (Runge-Kutta-Munthe-Kaas) methods which compose Runge-Kutta with the Lie-group exponential map.
- **Wendlandt, Marsden (1997)**, "Mechanical integrator derived from a discrete variational principle," *Physica D* **106**, 223–246. Earlier work establishing variational-integrator + discrete-variational-principle correspondence.
- **Kane, Marsden, Ortiz (1999)**, "Symplectic-energy-momentum preserving variational integrators," *J. Math. Phys.* **40**, 3353–3371. Relevant if the simultaneous-conservation question opens up in §5.

### §2.2 Lie-Poisson integrators (Touma-Wisdom / Breiter family)

**Foundational references:**
- **Ge, Marsden (1988)**, "Lie-Poisson Hamilton-Jacobi theory and Lie-Poisson integrators," *Phys. Lett. A* **133**, 134–139. The original construction.
- **Touma, Wisdom (1994)**, "Lie-Poisson integrator for rigid body dynamics in the solar system," *Astron. J.* **107**, 1189–1202. **Direct application to long-time solar-system rigid-body dynamics** — the closest astronomical-domain prior art outside the variational tier.
- **Breiter, Nesvorný, Vokrouhlický (2005)**, "Efficient Lie-Poisson integrator for secular spin dynamics of rigid bodies," *Astron. J.* **130**, 1267–1277. Modern refinement; explicitly designed for secular long-time integration of spin dynamics.

Lie-Poisson integrators preserve the **Casimir functions** of the Lie-Poisson bracket (for SO(3), this is the angular-momentum-squared invariant). They naturally preserve angular momentum to machine precision and conserve energy to bounded oscillation (backward error analysis). They do NOT necessarily preserve the symplectic 2-form, but the Casimir conservation is the closest analog for non-canonical Hamiltonian systems.

### §2.3 Splitting methods / non-canonical Hamiltonian systems (Hairer-Lubich-Wanner reference)

**Foundational reference:**
- **Hairer, Lubich, Wanner (2006)**, *Geometric Numerical Integration: Structure-Preserving Algorithms for Ordinary Differential Equations*, 2nd ed., Springer Series in Computational Mathematics 31. ISBN 9783642051579. **The authoritative textbook** on the whole field. Chapter IV covers "Conservation of First Integrals and Methods on Manifolds" (relevant to Lie-group methods); Chapter VII covers "Non-Canonical Hamiltonian Systems" (relevant to Lie-Poisson approaches); Chapter II/V cover symplectic/symmetric integrators in canonical Hamiltonian form. The book is cited explicitly as reference 7 in Lee-Leok-McClamroch (2007), confirming it is the standard reference both communities reach for.

**Splitting methods** decompose the Hamiltonian into two or more parts each individually integrable and compose the resulting exact flows. For the JWST-L2 problem, a natural split is **kinetic / gravitational** (the kinetic part is integrable on SO(3) × SO(3) per the Euler equations; the gravitational part is integrable in configuration space). The composition can be made symplectic to arbitrary order via Yoshida (1990) — *"Construction of high order symplectic integrators,"* Phys. Lett. A **150**, 262–268 — cited as reference 26 in Lee-Leok-McClamroch (2007). The 4th-order Yoshida split is a common practitioner choice.

**Splitting integrators for stochastic Lie-Poisson systems:** [arxiv:2111.07387](https://arxiv.org/pdf/2111.07387) (2022) explicitly demonstrates that splitting methods can preserve Casimirs in the stochastic-Lie-Poisson setting, suggesting the Phase 2 deterministic case has even cleaner Casimir-preservation properties.

### §2.4 Bou-Rabee-Marsden 2009 and stochastic variational integrators

**Bou-Rabee, Marsden (2009)**, "Hamilton-Pontryagin Integrators on Lie Groups," *Foundations of Computational Mathematics* **9**, 197–219. Generalises the Marsden-West variational-integrator framework to the Hamilton-Pontryagin variational principle (the constrained-Lagrangian form). Subsequent stochastic-variational follow-ups (Bou-Rabee, Owhadi 2010 and others) extend the framework to stochastic Hamiltonian systems for thermostat applications.

**Relevance to Phase 2:** **limited**. Phase 2 is deterministic — no SRP-as-stochastic-noise, no thermal-fluctuation forcing. The Bou-Rabee-Marsden 2009 deterministic extensions are subsumed by Lee-Leok-McClamroch (2007) for the JWST-L2 use case. Cite as a backstop in `refs.bib` but do not anchor on it.

### §2.5 Comparison table

| Property | Variational (LLM 2007) | Lie-Poisson (Touma-Wisdom 1994) | Splitting (Yoshida 1990 on HLW Ch VII) |
|---|---|---|---|
| Symplectic | Yes (discrete symplectic form) | Yes (Lie-Poisson sense; Casimir-preserving) | Yes (each step is a symplectic flow composition) |
| Momentum-preserving | Yes (discrete momentum map per Noether) | Yes (Casimirs are exactly preserved) | Yes (each split flow conserves its own integrals) |
| Energy behaviour at long times | Bounded oscillation (backward error analysis) | Bounded oscillation | Bounded oscillation |
| Group-structure preservation | Yes (without charts/reprojection) | Yes (Lie algebra dual; Casimirs) | Depends on each split-piece formulation |
| Convergence order | 1, 2 (typical); higher orders via composition | 2 (Touma-Wisdom standard); higher via composition | Even orders only via symmetric composition; 4th-order Yoshida common |
| Evaluations per step | 1 gravity/moment evaluation (LLM 2007 confirmed) | 1 (split into kinetic + force pieces) | n_split (number of composition stages) |
| Coupled SE(3) × SE(3) treatment | **Direct** — LLM 2007 paper covers this case explicitly | Solar-system spin dynamics; works for spin but treats translation separately | Natural: kinetic + gravity split; orbital + attitude coupling needs care |
| Closed-form Lagrange-top libration verification | Should match Ch 06 §6.4.3 to symplectic-bound accuracy | Should match; momentum-conserved | Should match; energy-bounded |
| Dzhanibekov regime numerical stability | Strong (no drift accumulation through saddle approach) | Strong (Casimir-conserved through saddle approach) | Strong but split-cadence-dependent |
| Implementation complexity for two-body coupled case | **Moderate-high** (LLM 2007 paper details the SE(3)×SE(3) implementation) | Moderate (Touma-Wisdom solar-system code as template) | Moderate (Hamiltonian split + Yoshida composition; well-documented in HLW Ch VII) |
| Direct prior-art applicability | **High** — LLM 2007 is JWST-L2 problem class | Moderate — Touma-Wisdom is single-body, asteroid-class | High — HLW Ch VII covers the theory; not the specific physics |

### §2.6 Recommendation

**Anchor Phase 2 on Lee-Leok-McClamroch (2007) variational integrator.** Three reasons:

1. **Direct prior art.** LLM 2007 applied this exact construction to the two-rigid-bodies-under-mutual-gravity problem. The dumbbell example in their paper is structurally identical to the JWST-L2 setup (two extended bodies with a coupled gravitational potential). No translation work is needed; the framework is on the shelf.
2. **Symplectic + momentum + group-structure preservation in one package.** The competing approaches preserve a subset (Lie-Poisson preserves Casimirs but not the symplectic form per se; pure-splitting preserves piece-wise integrals but the global preservation depends on composition design). Lee-Leok-McClamroch preserves all three structural properties through every step.
3. **Single gravity evaluation per step.** LLM 2007 explicitly states this is the efficiency win — the gravity computation is the expensive part in the JWST-L2 physics; saving evaluations matters.

**Secondary recommendation:** include a **splitting-method backstop** (e.g. Yoshida 4th-order composition of free-rotation + gravity-step) as a comparison integrator in `integrators/symplectic_se3_splitting.py`. The split formulation is simpler to author, easier to debug, and provides a sanity-check oracle for the variational integrator's outputs. If LLM 2007 implementation hits a snag, the splitting backstop becomes a viable Phase 2 ship.

**Decision boundary:** if the project owner declines Lee-Leok-McClamroch on grounds of implementation complexity, fall back to the 4th-order Yoshida splitting method. This is a recommended decision-point for §6 below.

---

## §3. Libration testbed design

### §3.1 Analytical oracle from Ch 06 §6.4.3

The chapter provides the closed-form Lagrange-top + gravity-gradient libration solution via Podolský 2018 §8.3 effective-potential method, with substitution `g → 3Gm_B / |ρ|^3`. The key analytical predictions per Ch 06 §6.4.3:

- **Boxed cyclic-coordinate first integrals** L_ψ, L_φ (Podolský Eqs. 8.10–8.11 analogs)
- **Boxed effective potential** V_ef(ϑ)
- **Boxed quadrature for ϑ(t)** 
- **Three precession regimes** (a)/(b)/(c) per Lagrange-top classical taxonomy
- **Finite-amplitude period** T_lib(ϑ_amp) = (2/π) K(k) · T_lib^linear, where K(k) is the complete elliptic integral of the first kind with modulus k determined by the libration amplitude

The testbed validates the Phase 2 integrator by:
1. Initialising in the small-amplitude libration regime (where the linear period prediction holds within ~1 %)
2. Integrating for ~100 libration periods (to test the symplectic-vs-RK4 difference becomes visible; RK4 would drift the libration phase by O(0.01 rad) over 100 periods at typical Phase 1 timesteps, symplectic should hold to within the rounding-error budget)
3. Extracting the integrated libration period (FFT of the libration coordinate trace, e.g. `θ - θ_eq`)
4. Comparing against the analytical T_lib prediction

**Pass gate:** integrated period within 0.1 % of analytical for 100 periods at default timestep. Cross-check at finer timesteps (Δt/2, Δt/4) confirms order-of-convergence is consistent with the integrator's claimed order (e.g. ~2 for LLM 2007 standard integrator, ~4 for Yoshida-composed variant).

### §3.2 JWST-specific libration considerations

JWST's actual attitude dynamics at L2 are dominated by **solar-radiation-pressure torque on the sunshield**, not by **gravity-gradient torque** as the prompt-mode coupling. Per [NASA NTRS 20140008865 (Yoon, Rosales, Richon 2014)](https://ntrs.nasa.gov/citations/20140008865), the "major challenge to accurate JWST OD during the nominal science phase results from the unusually large solar radiation pressure force acting on the huge sunshield." Per [JWST docs (STSCI)](https://jwst-docs.stsci.edu/jwst-observatory-characteristics-and-performance/jwst-orbit), "the solar torque is balanced by reaction wheels, but periodically, the accumulated momentum is dumped by firing thrusters" — confirming SRP is the dominant attitude-disturbance source at L2.

**Implication for Phase 2:** the libration testbed should use a **JWST-like-but-idealised** configuration — geometry matching `make_jwst_like()` per Phase 1, but with SRP forcing **turned off** (deferred to Phase 3 per the existing SONNET-BRIEF.md §1). This isolates the gravity-gradient libration dynamics that the analytical oracle covers, without polluting the validation with SRP physics the analytical oracle doesn't model.

Verify the gravity-gradient libration period for the JWST-like reference body at the Phase 1 default ρ (Sun-Earth L2 separation, 1.5 × 10⁹ m per STSCI docs) gives a libration period in the band where the testbed is observable — the §6.4.3 numerical sanity check predicts ~66 day libration period for the JWST-like prolate composite at L2 against the gravity-gradient torque alone, well within the testbed's 100-period observation window.

### §3.3 Alternative gravity-gradient libration anchors

The classical spacecraft-attitude-dynamics literature provides additional libration anchors worth checking against:

- **Hughes (1986)**, *Spacecraft Attitude Dynamics*, John Wiley & Sons (Dover reprint 2004). Chapters 7–8 on gravity-gradient libration with rate-and-phase decompositions. Practitioner standard.
- **Wertz (1978, ed.)**, *Spacecraft Attitude Determination and Control*, Kluwer. Reference handbook; cross-check the §3.1 analytical predictions against Wertz's Eq. 18.13 (planet-pointed libration period) for sanity.
- **Curtis (2014)**, *Orbital Mechanics for Engineering Students*, 3rd ed., Elsevier. Pedagogical exposition of gravity-gradient stabilisation.

These are **not** Phase 2 critical-path; cite in `refs.bib` and use as cross-references if the §3.1 testbed gives unexpected results.

---

## §4. Lyapunov-spectrum methodology

### §4.1 The Benettin-Galgani-Giorgilli-Strelcyn 1980 algorithm

**Foundational reference:**
- **Benettin, Galgani, Giorgilli, Strelcyn (1980)**, "Lyapunov characteristic exponents for smooth dynamical systems and for Hamiltonian systems; a method for computing all of them. Part 1: Theory; Part 2: Numerical application," *Meccanica* **15**, 9–20 and 21–30. The original two-part paper establishing the still-standard reorthogonalisation method.

The algorithm: propagate the full state vector AND a basis of tangent-space perturbation vectors. After every fixed-interval τ_QR, Gram-Schmidt-orthonormalize the perturbation basis (equivalently QR-decompose), accumulate the logarithms of the diagonal entries of R, and divide by total elapsed time to estimate the Lyapunov exponents.

Per the search results (Skokos 2010 chapter, Geist-Parlitz-Lauterborn 1990 review):
- **Gram-Schmidt** is the original method (Shimada-Nagashima, Benettin et al.)
- **Householder transformations** suggested by Eckmann-Ruelle as numerically superior
- **QR or SVD decomposition** are the modern continuous-time variants
- Multiple modern algorithms (Wolf 1985, Geist et al. 1990, Skokos 2010, the differential algorithm of arxiv:1008.3368 2010, and the "fastest simplified method" of Springer Nonlinear Dyn 2018) all start from the Benettin et al. backbone

**Standard reference for review:** [Skokos C. (2010), "The Lyapunov Characteristic Exponents and Their Computation," in *Dynamics of Small Solar System Bodies and Exoplanets*, Lect. Notes Phys. 790, Springer.](https://link.springer.com/chapter/10.1007/978-3-642-04458-8_2) Authoritative modern review with all the algorithm details.

### §4.2 Implementation specifics for SO(3) × SO(3)

For Phase 2, the tangent-space cocycle propagation must respect the Lie-group structure. Two approaches:

1. **Variational equations in Euclidean form:** linearise the equations of motion in the local chart (Euler-angle or quaternion-component space), propagate the tangent vectors in R^n, reorthogonalize in R^n. Simpler to author; loses some Lie-group geometric properties at the cocycle step.

2. **Lie-algebraic variational propagation:** propagate the tangent vectors in the Lie algebra `so(3) × so(3) ≅ R^6`, applying the body-frame variational equations directly. Reorthogonalisation happens in R^6 (the Lie algebra is a vector space, so Gram-Schmidt works natively). More geometrically clean; pairs naturally with the LLM 2007 variational integrator framework.

**Recommendation:** start with approach 1 (Euclidean variational equations in quaternion-component coordinates) as a quick first-cut implementation. Migrate to approach 2 only if a numerical-stability issue surfaces (e.g. the tangent vectors drift off the Lie algebra unphysically).

### §4.3 QR cadence

The QR-decomposition cadence τ_QR is the central tunable parameter:
- **Too small** (e.g. τ_QR < timestep): wasted compute on re-orthogonalisation that produces near-zero norm change
- **Too large** (e.g. τ_QR > 1/λ_max, where λ_max is the largest exponent): the perturbation basis becomes numerically singular (all vectors align with the most-expanding direction) before reorthogonalisation; QR factorisation loses precision

**Standard guideline** (per Skokos 2010 and Geist-Parlitz-Lauterborn 1990): τ_QR ≈ 0.1 / λ_max_expected. For the Dzhanibekov regime, λ_max scales with the asymmetry parameter; order-of-magnitude estimate 1e-3 to 1e-1 s^-1 for moderate-asymmetry bodies. So τ_QR in the range 1–100 s is a reasonable starting point.

**Empirical autocalibration:** run a short pilot integration, estimate λ_max from the unnormalised growth rate, set τ_QR = 0.1 / λ_max_pilot, re-run for production. Document in the testbed report.

### §4.4 Finite-time vs asymptotic estimators

Per Geist-Parlitz-Lauterborn (1990): the Lyapunov-exponent estimate from a finite-time integration converges to the asymptotic value at rate ~1/√T (for chaotic systems with smooth invariant measure). To get exponent estimates to within 1% relative error, you need T such that √T × λ ~ 100 — i.e. integration time T ~ 10^4 / λ². For λ ~ 0.01 s^-1, that's T ~ 10^8 s (~3 years). **This is the canonical "Lyapunov-spectrum needs long integration" budget.**

Phase 2 should NOT aim for asymptotic convergence on a Sonnet-session timescale. Instead, the Phase 2 acceptance gate should be:
- **Qualitative**: λ_max > 0 in the Dzhanibekov regime; λ_max ≈ 0 in the symmetric-top regime; ordering of exponents respects the Hamiltonian-system constraint (each positive λ_i paired with a negative −λ_i, plus two zero exponents from the Hamiltonian + reduced-Hamiltonian invariants)
- **Order-of-magnitude quantitative**: λ_max matches the analytical Dzhanibekov-flip-frequency-derived expectation per Ch 06 §6.5.3 + Landau-Lifshitz §37 closed-form to within a factor of ~3
- **Symmetric-top quantitative**: λ_max ≈ 0 to machine precision (because the symmetric top is integrable, the Lyapunov spectrum is all-zero in exact arithmetic)

### §4.5 Dzhanibekov-regime Lyapunov signature

The Dzhanibekov / tennis-racket / intermediate-axis instability is **structurally** chaotic in the sense that small perturbations from the unstable equilibrium grow exponentially. But the dynamics also have a **conserved quantity** (the total angular momentum magnitude L = |L|), which means the trajectory is confined to a 2-sphere in body-frame angular-velocity space, and the dynamics on the sphere are integrable (Euler-equation solution in terms of Jacobi-elliptic functions per Landau-Lifshitz §37).

**This integrability is the key subtlety.** Per [arxiv:2003.13539 (2020)](https://arxiv.org/abs/2003.13539) "Geometric Origin of the Tennis Racket Effect":
> "The tennis racket effect is a geometric phenomenon which occurs in a free rotation of a three-dimensional rigid body. In a complex phase space, we show that this effect originates from a pole of a Riemann surface and can be viewed as a result of the Picard-Lefschetz formula. We prove that a perfect twist of the racket is achieved in the limit of an ideal asymmetric object. We give upper and lower bounds to the twist defect for any rigid body, which reveals the robustness of the effect."

So the **free-rotation Dzhanibekov effect is integrable** (Jacobi-elliptic + Riemann-surface poles). It is **not chaotic** in the Lyapunov-positive-exponent sense in vacuum. **Real chaos** appears only when the system is **perturbed off the integrable manifold** — by gravity-gradient coupling to the orbital state, by SRP forcing, by mass asymmetry breaking the angular-momentum conservation, etc.

**For the Phase 2 Lyapunov testbed, this means:**
- The pure torque-free asymmetric top has λ_max = 0 (integrable). The §6.5.3 Jacobi-elliptic parametrisation is the analytical oracle.
- The gravity-gradient-coupled JWST-like body in orbit around the gravitational primary breaks integrability via the orbital-attitude coupling. λ_max should be > 0 — but possibly small. This is the **physically interesting Lyapunov regime** for the JWST-L2 problem.
- A useful **calibration testbed** would integrate the same asymmetric body under three scenarios: (a) free in vacuum (λ_max should be ≈ 0 per Jacobi-elliptic integrability); (b) in gravity-gradient coupling at a stationary equilibrium (λ_max should be small per libration-mode-only dynamics); (c) in gravity-gradient coupling near a Dzhanibekov-flip initial condition (λ_max should be > 0 and reflect the orbital-attitude coupling-induced sensitivity).

### §4.6 Recent dynamical-systems treatments to cite

- [arxiv:2003.13539](https://arxiv.org/abs/2003.13539) — Geometric origin via Picard-Lefschetz; explicit twist-defect bounds for any rigid body
- [arxiv:1606.08237](https://arxiv.org/abs/1606.08237) — "The tennis racket effect in a three-dimensional rigid body" — earlier paper in the same line (fetch attempted but timed out; cite by URL only)
- **Murakami et al.** — various papers in *Physica D* on the tennis-racket Hamiltonian structure
- **Trivailo, Kojima** — recent work on the Dzhanibekov dynamics with applications to spacecraft attitude control

These are NOT Phase 2 implementation references; they motivate the §4.5 framing of why the Lyapunov spectrum in the free-rotation case is zero, and what the JWST-coupled case adds.

---

## §5. C2 cross-test tolerance calibration for symplectic integrators

### §5.1 The Phase 1 tolerance

Phase 1's C2 cross-tests used **Higham §4.2 Option G** ULP×√N bounds. The §4.2 Eq. 4.6 form is:

  | δx | ≤ γ_n · |x|

where γ_n = nu / (1 − nu), u is the unit roundoff (~1.1e-16 for double), and n is the operation count. For ~10^4 operations per timestep × 10^4 timesteps = 10^8 ops total, γ_n is bounded by ~10^8 × u ~ 10^-8 in worst case. For Option G specifically — the probabilistic version with √N replacing N — the worst-case constant becomes √N × u ~ 10^4 × 1.1e-16 ~ 10^-12.

Phase 1 trajectory cross-tests landed at worst-rel ~3.1e-11 — about an order of magnitude above the Option G bound, well within the gate ceiling.

### §5.2 What changes for symplectic integrators

Phase 1's RK4 has **discretisation error per step ~Δt^4** but **no qualitative structure preservation**: energy and momentum both drift. The trajectory cross-test thus directly compares numerical output to a high-resolution reference trajectory — the comparison is *trajectory-vs-trajectory*, and both numerical and reference share the same Δt^4 discretisation error scaling.

For Phase 2's symplectic integrator, the situation is different:
- **Energy** is bounded in oscillation (backward error analysis: the symplectic integrator solves an exact Hamiltonian system whose Hamiltonian differs from the original by O(Δt^p), where p is the integrator order — Hairer-Lubich-Wanner 2006 Ch IX). So energy oscillates around the true value with bounded amplitude.
- **Momentum** is exactly conserved (per the discrete-Noether theorem — Marsden-West 2001 §4).
- **Trajectory** still has O(Δt^p) discretisation error per step that accumulates as Δt^p × T_total.

So a trajectory C2 cross-test of two implementations of the same symplectic integrator (e.g., the Python ref-impl vs a future Fortran port) at the same step-size should give:
- Energy conservation: identical to ULP precision (both compute the same quantity)
- Momentum conservation: identical to ULP precision
- Trajectory: should match to ULP precision if both implementations are bit-identical floating-point operations (since the integrator step is a deterministic floating-point function)

**Recommended tolerance for Phase 2 symplectic C2 cross-tests:**

| Quantity | Phase 1 (RK4) tolerance | Phase 2 (symplectic) tolerance |
|---|---|---|
| Energy at integration endpoint | ULP×√N (Option G) | ULP×√N (Option G); same |
| Angular momentum magnitude at endpoint | ULP×√N | ULP precision (because momentum is exactly conserved in the integrator) |
| Trajectory state vector at endpoint | ULP×√N | ULP×√N — same scaling, because the per-step round-off model is the same regardless of integrator structure |

**Net change:** angular-momentum cross-test tightens from `ULP×√N` to **ULP precision** (no √N safety factor). Energy and trajectory tolerances remain the same Higham §4.2 Option G form.

### §5.3 Long-time-window cross-tests for Phase 2

Phase 1 trajectory cross-tests run for ~10^4 steps. Phase 2 should add a **long-time-window cross-test** at ~10^6 steps to validate the symplectic property:

- **Phase 1 RK4 would fail** this test (energy drifts by O(Δt^4 × 10^6 × T_step) > Option G bound at 10^6 steps).
- **Phase 2 symplectic should pass** this test (energy oscillates within bounded amplitude regardless of integration length).

This is itself a **conservation-property test**, not a trajectory-cross-test. The acceptance criterion is qualitative ("energy excursion bounded vs time" plot is flat-mean-with-bounded-oscillation rather than drifting upward) plus a quantitative bound on the energy excursion amplitude (per Hairer-Lubich-Wanner backward-error-analysis: amplitude ~Δt^p × ω where ω is a characteristic system frequency).

### §5.4 Quantitative backward-error-analysis amplitude bound

Per Hairer-Lubich-Wanner 2006 Ch IX Theorem 8.1 (approximate; the exact statement requires more care): a symplectic integrator of order p applied to a Hamiltonian system H produces a numerical solution that exactly satisfies a modified Hamiltonian H̃ with H̃ − H = O(Δt^p × C(H)) for some constant C depending on the smoothness of H.

For Phase 2's symplectic integrator at order p (typically 2 for the basic Lee-Leok-McClamroch, 4 for the Yoshida-composed variant), and for the JWST-L2 Hamiltonian, the predicted energy-excursion amplitude is:

  ΔH_amplitude / H ~ (Δt × ω_max)^p

where ω_max is the largest characteristic frequency of the system (for JWST-L2, the orbital angular velocity at L2 is ω_orbital ≈ 2π/(6 months) ≈ 4e-7 rad/s, but the rigid-body rotation at JWST design spin rate is much faster — typically 1e-4 to 1e-2 rad/s depending on whether the body is on the spin-stabilised vs three-axis-stabilised side). For a representative Δt = 1 s and ω_max = 1e-2 rad/s, the predicted amplitude is ΔH/H ~ (1e-2)^2 = 1e-4 for order 2, or 1e-8 for order 4.

**Use this as the Phase 2 acceptance criterion** for the long-time-window energy-conservation test: energy excursion stays bounded within the predicted amplitude over ~10^6 steps, with no secular drift visible on a residual-vs-time plot.

---

## §6. Open questions + decision points for project owner

The following points need project-owner direction before the Sonnet executor fires the brief:

**OQ-PHASE2-1: Variational vs Yoshida-splitting integrator choice.** §2.6 recommends anchoring on Lee-Leok-McClamroch (2007) variational integrator (with Yoshida-splitting as backstop). **Decision:** confirm LLM 2007 as the primary integrator, OR fall back to Yoshida-splitting as primary, OR implement both with the LLM as primary.

**OQ-PHASE2-2: Lie-algebraic vs Euclidean tangent-space Lyapunov propagation.** §4.2 recommends starting with Euclidean (quaternion-component) variational equations + Gram-Schmidt and migrating to Lie-algebraic only if needed. **Decision:** confirm Euclidean-first OR mandate Lie-algebraic from the start.

**OQ-PHASE2-3: Lyapunov long-time-window budget.** §4.4 notes the canonical Lyapunov-spectrum convergence budget is T ~ 10^4/λ² seconds. The Sonnet executor needs an upper bound on integration time per Lyapunov-spectrum-extraction run. **Decision:** specify maximum wall-time per integration run (e.g. 1 h, 4 h, 8 h) so the Sonnet executor can size Δt and the number of integration windows accordingly.

**OQ-PHASE2-4: SRP forcing scope.** §3.2 confirms SRP is JWST's dominant attitude-disturbance source. Phase 2's existing SONNET-BRIEF.md defers SRP to Phase 3. **Decision:** confirm Phase 2 has SRP forcing OFF; OR pull SRP into Phase 2 scope (would significantly expand the brief).

**OQ-PHASE2-5: Fortran port timing.** Phase 2 SONNET-BRIEF.md §1 currently states "Phase 2 does NOT touch the Fortran backend at this milestone." The Python-only Phase 2 + Fortran-port-in-v0.2.x cadence matches the Phase 1 pattern. **Decision:** confirm Python-only Phase 2 (Fortran port queues for v0.2.1+), OR include Fortran port in Phase 2 scope.

**OQ-PHASE2-6: C2 cross-test tolerance recalibration.** §5.2 recommends tightening the angular-momentum C2 tolerance from `ULP×√N` to `ULP precision` for the symplectic integrator. **Decision:** confirm tolerance tightening, OR keep the `ULP×√N` Phase 1 convention.

**OQ-PHASE2-7: Bibliography additions to refs.bib.** This landscape report cites ~15 new references not in Phase 1's `refs.bib`. **Decision:** confirm which references get added to refs.bib; the Sonnet executor will update accordingly.

---

## §7. Recommended Sonnet brief amendments

The existing SONNET-BRIEF.md at `_handoffs/phase-2-symplectic-lie-group/SONNET-BRIEF.md` is already a serviceable implementation contract. Based on this landscape's findings, recommended amendments:

### §7.1 Add a §0.5 "research-context anchor" pointing at this landscape

One-paragraph block near the brief's top:

> **Research context.** Phase 2 algorithmic + literature triangulation lives in [`../../_specs/PHASE-2-INTEGRATOR-LANDSCAPE-v0.1.md`](../../_specs/PHASE-2-INTEGRATOR-LANDSCAPE-v0.1.md). Read it first before opening the brief; the integrator-choice rationale, the libration-testbed analytical-oracle calibration, the Lyapunov-methodology specifics, and the C2 tolerance recalibration are all derived there.

### §7.2 Update §1 (Scope summary) integrator-choice cell

Currently: "RKMK / discrete-Euler-Poincaré integrator on SO(3) × SO(3)". Per §2.6 recommendation:

> **Symplectic Lie-group variational integrator** — anchor on Lee-Leok-McClamroch (2007) construction for the full body problem in orbital mechanics, applied to SO(3) × SO(3) attitude dynamics with two-body translation handled separately. Implementation backstop: 4th-order Yoshida-composition splitting of kinetic + gravitational pieces per Hairer-Lubich-Wanner 2006 Ch VII. **Decision pending on OQ-PHASE2-1.**

### §7.3 Tighten §2.1 module layout

Currently includes `symplectic_so3.py` and `symplectic_se3.py`. Per §2.6 backstop recommendation, add:

```
integrators/
├── __init__.py
├── rk4.py                              # extracted from current dynamics.py; preserved
├── symplectic_se3_variational.py       # NEW — Lee-Leok-McClamroch 2007 anchor
└── symplectic_se3_splitting.py         # NEW — Yoshida 1990 backstop
```

### §7.4 Add §3.1bis (libration-testbed analytical-oracle calibration)

Reference Ch 06 §6.4.3 closed-form Lagrange-top libration explicitly. Pass-gate: integrated period within 0.1 % of analytical over 100 libration periods. See landscape §3.1.

### §7.5 Add §4.5 Lyapunov-spectrum testbed clarification

The torque-free asymmetric top is integrable (λ_max = 0); the gravity-gradient-coupled case is what generates Lyapunov-positive chaos. See landscape §4.5 for the calibration scenarios (free vs gravity-gradient-stationary vs gravity-gradient-Dzhanibekov-IC).

### §7.6 Add §5.2bis C2 tolerance recalibration table

Per landscape §5.2 — angular momentum tightens to ULP precision; energy + trajectory stay at Option G. Decision pending on OQ-PHASE2-6.

### §7.7 Add §6.1 OQ-PHASE2-1..7 to the brief's open-questions section

Lift the §6 OQ list from this landscape into the brief verbatim. Resolve before Sonnet fires.

### §7.8 Recommended brief acceptance gates (verbatim, for §10 of SONNET-BRIEF.md)

| Gate | Criterion |
|---|---|
| **AG-INT-1** | Lee-Leok-McClamroch (2007) variational integrator implemented per §2.6; Yoshida-splitting backstop implemented per §2.6 if OQ-PHASE2-1 directs both |
| **AG-INT-2** | Energy excursion over 10^6-step integration bounded within predicted backward-error-analysis amplitude per landscape §5.4 |
| **AG-INT-3** | Angular momentum conserved to ULP precision at integration endpoint per landscape §5.2 |
| **AG-LIB-1** | Libration testbed matches Ch 06 §6.4.3 analytical period within 0.1 % over 100 libration periods per landscape §3.1 |
| **AG-LIB-2** | Libration-period order-of-convergence test confirms integrator's claimed order at Δt/2 and Δt/4 |
| **AG-LYAP-1** | Torque-free asymmetric-top Lyapunov spectrum has λ_max ≈ 0 to machine precision (integrability check per landscape §4.5) |
| **AG-LYAP-2** | Gravity-gradient-Dzhanibekov-IC Lyapunov spectrum has λ_max > 0 and matches Ch 06 §6.5.3 + Landau-Lifshitz §37 closed-form characteristic frequency to within factor of ~3 |
| **AG-LYAP-3** | Hamiltonian-system Lyapunov-spectrum pairing constraint (each +λ paired with −λ, plus zero modes) is observable in the computed spectrum |
| **AG-NUM** | All numerical claims in the chapter-equivalent docs (canonical or release-notes) are grep-traceable to source script outputs |
| **AG-VOICE** | If Phase 2 outputs chapter-equivalent docs (canonical-tier Ch 07 or 08), they pass canonical-tier voice convention per project standard (not Shad-tier) |
| **AG-SANITISE** | Sanitisation pre-flight CLEAN per `_config/SANITISATION-POLICY-v0.1.md` |
| **AG-BIB** | refs.bib updated with Lee-Leok-McClamroch 2007, Marsden-Pekarsky-Shkoller 1999/2000, Marsden-West 2001 Acta Numerica, Iserles-Munthe-Kaas-Nørsett-Zanna 2000 Acta Numerica, Benettin-Galgani-Giorgilli-Strelcyn 1980, Skokos 2010, Touma-Wisdom 1994, Yoshida 1990, Bou-Rabee-Marsden 2009 per OQ-PHASE2-7 |

---

## §8. References (consolidated)

Primary references for direct anchoring:

- Lee, T., Leok, M., McClamroch, N.H. (2007), "Lie group variational integrators for the full body problem in orbital mechanics," *Celestial Mechanics and Dynamical Astronomy* **98**, 121–144. DOI: [10.1007/s10569-007-9073-x](https://doi.org/10.1007/s10569-007-9073-x). ([SpringerLink](https://link.springer.com/article/10.1007/s10569-007-9073-x))
- Marsden, J.E., Pekarsky, S., Shkoller, S. (1999), "Discrete Euler-Poincaré and Lie-Poisson equations," *Nonlinearity* **12**, 1647–1662.
- Marsden, J.E., Pekarsky, S., Shkoller, S. (2000), "Symmetry reduction of discrete Lagrangian mechanics on Lie groups," *J. Geom. Phys.* **36**, 139–150.
- Marsden, J.E., West, M. (2001), "Discrete mechanics and variational integrators," *Acta Numerica* **10**, 357–514.
- Hairer, E., Lubich, C., Wanner, G. (2006), *Geometric Numerical Integration: Structure-Preserving Algorithms for ODEs*, 2nd ed., Springer Series in Computational Mathematics 31. ISBN 9783642051579.
- Iserles, A., Munthe-Kaas, H.Z., Nørsett, S.P., Zanna, A. (2000), "Lie-group methods," *Acta Numerica* **9**, 215–365.
- Benettin, G., Galgani, L., Giorgilli, A., Strelcyn, J.-M. (1980), "Lyapunov characteristic exponents for smooth dynamical systems and for Hamiltonian systems; a method for computing all of them. Part 1: Theory; Part 2: Numerical application," *Meccanica* **15**, 9–20 and 21–30.
- Skokos, C. (2010), "The Lyapunov Characteristic Exponents and Their Computation," in *Dynamics of Small Solar System Bodies and Exoplanets*, Lect. Notes Phys. **790**, Springer. ([SpringerLink](https://link.springer.com/chapter/10.1007/978-3-642-04458-8_2))
- Touma, J., Wisdom, J. (1994), "Lie-Poisson integrator for rigid body dynamics in the solar system," *Astron. J.* **107**, 1189–1202.
- Yoshida, H. (1990), "Construction of high order symplectic integrators," *Phys. Lett. A* **150**, 262–268.
- Bou-Rabee, N., Marsden, J.E. (2009), "Hamilton-Pontryagin integrators on Lie groups," *Foundations of Computational Mathematics* **9**, 197–219.

Secondary references for context:

- Lee, T., Leok, M., McClamroch, N.H. (2005), "A Lie group variational integrator for the attitude dynamics of a rigid body with applications to the 3D pendulum," *Proc. IEEE Conf. Control Applications*, pp. 962–967. ([PDF](https://mathweb.ucsd.edu/~mleok/pdf/LeLeMc2005_lgvi.pdf))
- Lee, T., Leok, M., McClamroch, N.H. (2007), "Lie group variational integrators for the full body problem," *Comput. Methods App. Mech. Eng.* (companion paper to the CMDA one).
- Wendlandt, J.M., Marsden, J.E. (1997), "Mechanical integrator derived from a discrete variational principle," *Physica D* **106**, 223–246.
- Kane, C., Marsden, J.E., Ortiz, M. (1999), "Symplectic-energy-momentum preserving variational integrators," *J. Math. Phys.* **40**, 3353–3371.
- Ge, Z., Marsden, J.E. (1988), "Lie-Poisson Hamilton-Jacobi theory and Lie-Poisson integrators," *Phys. Lett. A* **133**, 134–139.
- Breiter, S., Nesvorný, D., Vokrouhlický, D. (2005), "Efficient Lie-Poisson integrator for secular spin dynamics of rigid bodies," *Astron. J.* **130**, 1267–1277.
- Scheeres, D.J. (2002), "Stability in the Full Two-Body Problem," *Celestial Mech. Dyn. Astr.* **83**, 155–169.
- Geist, K., Parlitz, U., Lauterborn, W. (1990), "Comparison of Different Methods for Computing Lyapunov Exponents," *Progress of Theoretical Physics* **83**, 875–893. ([PDF](http://carretero.sdsu.edu/teaching/M-638/lectures/ProgTheorPhys_1990-Geist-875-93.pdf))

Recent dynamical-systems treatments of Dzhanibekov / tennis-racket effect:

- "Geometric Origin of the Tennis Racket Effect" (2020) — arxiv:[2003.13539](https://arxiv.org/abs/2003.13539)
- "The tennis racket effect in a three-dimensional rigid body" — arxiv:[1606.08237](https://arxiv.org/abs/1606.08237)

JWST orbital + attitude context:

- Yoon, S., Rosales, J., Richon, K. (2014), "James Webb Space Telescope Orbit Determination Analysis," NASA NTRS [20140008865](https://ntrs.nasa.gov/citations/20140008865), International Symposium on Space Flight Dynamics, Laurel MD.
- STSCI, "JWST Orbit," *JWST User Documentation*, [jwst-docs.stsci.edu/jwst-observatory-characteristics-and-performance/jwst-orbit](https://jwst-docs.stsci.edu/jwst-observatory-characteristics-and-performance/jwst-orbit) (accessed 2026-06-01).

Spacecraft attitude dynamics classical references:

- Hughes, P.C. (1986), *Spacecraft Attitude Dynamics*, John Wiley & Sons (Dover reprint 2004).
- Wertz, J.R. (ed.) (1978), *Spacecraft Attitude Determination and Control*, Kluwer.
- Curtis, H. (2014), *Orbital Mechanics for Engineering Students*, 3rd ed., Elsevier.

Higham tolerance baseline:

- Higham, N.J. (2002), *Accuracy and Stability of Numerical Algorithms*, 2nd ed., SIAM. §4.2 Eq. 4.6 + Option G probabilistic ULP×√N bound — the Phase 1 + Phase 2 C2 cross-test tolerance baseline.

Recent splitting + Lie-Poisson literature:

- "Splitting integrators for stochastic Lie-Poisson systems" (2022) — arxiv:[2111.07387](https://arxiv.org/abs/2111.07387)
- "A differential algorithm for the Lyapunov spectrum" (2010) — arxiv:[1008.3368](https://arxiv.org/abs/1008.3368)

---

**End PHASE-2-INTEGRATOR-LANDSCAPE-v0.1.md.**  
*Landscape v0.1 informs the Phase 2 Sonnet brief; project owner resolves §6 OQ-PHASE2-1..7 before brief fires.*
