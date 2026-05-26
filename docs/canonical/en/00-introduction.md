# lege-artis/jwst-l2-coupled-dynamics — canonical tier — Introduction (EN)

> **Audience.** Mathematician, physicist, academic reviewer, contributor wanting to understand *why* the coupled-dynamics code in this project computes what it claims to compute, derived from first principles in modern frame-invariant form.
>
> **Reading time.** ~30 minutes for this Introduction. ~3-4 hours for the full canonical tier (this + chapters 01–06) once complete.
>
> **License.** CC-BY-SA-4.0 (intended; project not yet publicly mirrored).
>
> **Companion tiers.** Engineer-tier `docs/engineer/en/00-quick-start.md` (forthcoming) will cover the same scope from the practitioner-using-the-library angle. Shad-tier `SHAD-NARRATIVE.md` (already shipped) covers the engineer-narrated walkthrough of the first-cut prototype.

---

## §1. Scope of the canonical tier

This project — currently scaffolded under `experiments/jwst-l2-first-cut/` in the workspace, intended as the seed of future `lege-artis/jwst-l2-coupled-dynamics` — provides a reference implementation of two-rigid-body gravitationally-coupled dynamics, anchored on a JWST-like primary and a probe-like secondary in the vicinity of (eventually) the Earth–Sun L2 equilibrium.

The canonical tier derives, from first principles, the equations of motion that govern this system. The derivation proceeds in modern frame-invariant language and specialises to specific frames only when calculation demands it.

The **scope of v0.1 of this canonical tier** (the scope authored here) is:

| Domain element | Treatment |
|---|---|
| Two rigid bodies (extended, with full inertia tensors) | Yes |
| Mutual gravitational potential between them | Yes, derived from first principles (multipole expansion) |
| Rigid rotation about the centre of mass | Yes, derived via Lagrange + Euler-Poincaré reduction |
| Relative orbital motion in the centre-of-mass frame | Yes |
| Restricted three-body L2 dynamics (Sun-Earth-spacecraft) | **Deferred to v0.2** of this canonical tier |
| Relativistic 1PN/2PN corrections | **Explicitly excluded** (justified in §6 below; relevant velocity ratio is $v/c \sim 10^{-5}$ for JWST orbital scales, contribution to rotational coupling negligible) |
| External perturbations (solar radiation pressure, third-body Moon perturbation, etc.) | **Out of scope** at canonical-tier; engineer-tier will catalogue them |
| Internal flexibility / vibration modes | **Out of scope**; rigid-body assumption is the canonical-tier baseline |

The **v0.1 chapter set** (the chapters that complete this tier at v0.1):

| Chapter | Title | Status | CS-language anchor (Podolský 2018) |
|---|---|---|---|
| 00 | Introduction (this file) | **Authored** | — (introductory; CS translation kicks off in `_handoffs/cs-canonical-tier-v0.1/`) |
| 01 | Configuration manifold and variational principle | **Scaffolded** | §6.2 + §6.4 (orthogonal-matrix + Euler-angle chart as elementary precursor; Lie-group abstraction is modern Marsden-Ratiu, not in Podolský) |
| 02 | Kinetic energy on $SE(3) \times SE(3)$ in intrinsic form | **Authored** (v0.1, round 2 augmentations 2026-05-24 AM: §2.5.1 d**I**/dt sensitivity, §2.6.1 CoM symplectic, §2.7.1 left-invariance) | §6.3 Eq. (6.4) body-frame Ω + (6.8) boxed transport equation + (6.2)/(6.10) trivialisation; §7.1 (7.3)/(7.7)/(7.9)/(7.11) inertia tensor + (7.12) Steinerova věta |
| 03 | Mutual gravitational potential and multipole expansion | **Authored** (v0.1, with §3.7.2 OQ-FORT-1 sign-fix landed 2026-05-24 PM) | Partial: §7.1 inertia-tensor algebra as bridge to gravitational quadrupole (§3.5); §8.3 Lagrange-top libration as classical precedent for gravity-gradient libration (§3.7). Modern multipole + gravity-gradient material is outside Podolský's scope. |
| 04 | Euler–Lagrange equations and Euler–Poincaré reduction | **Authored** (v0.1, canonical-tier round 3 2026-05-25: §4.2 Lagrangian assembly + §4.3 Approach A Euler-Lagrange in Euler-angle chart + §4.4 Approach B Euler-Poincaré reduction with explicit $\mathrm{ad}^*$ derivation + §4.5 cross-validation A≡B + §4.6 translational equation + back-reaction force + §4.7 boxed total equation set + line-by-line cross-validation against `dynamics.py`; all four chapter-internal OQs 4.1-4.4 resolved; forward-looking OQ-4.5..4.8 surfaced for v0.1.1/v0.2) | §7.2 boxed (7.15) Euler dynamic equations + §7.3 full Lagrangian-formalism derivation via Euler angles + §6.4 (6.14)–(6.17) chart machinery — direct CS anchor for Approach A; the two-in-parallel methodology this chapter formalises is implicit in Podolský's §7.2 + §7.3 two-derivation presentation. Euler-Poincaré reduction (Approach B) is modern Marsden-Ratiu, not in Podolský. |
| 05 | Body-frame Euler equations as a specialisation | **Authored** (v0.1, canonical-tier round 3 2026-05-25: §5.2 principal-axis specialisation with component-by-component expansion + §5.3 symmetry-class worked paragraphs + §5.4 gravity-gradient torque substituted in principal-axis form + §5.4.1 libration-equilibrium cross-check vs Ch 03 §3.7.2 + §5.5 conservation laws torque-free and torqued + §5.6 testbed pointers; all four chapter-internal OQs 5.1-5.4 resolved; forward-looking OQ-5.5..5.6 for v0.2) | §7.2 boxed (7.15) — direct CS analog of the chapter's central result; §8.1 (8.1)–(8.2) free symmetric top as analytical testbed; §8.2 setrvačník se třením as damped-attitude variant |
| 06 | Analytical limits and validation gates | **Authored** (v0.1, canonical-tier round 3 2026-05-25: §6.1 Class A/B discipline + §6.2 Kepler-limit testbed + §6.3 free symmetric top testbed (implemented in first-cut) + §6.4 gravity-gradient libration testbed (Class A linearised + Class B nonlinear) + §6.5 Dzhanibekov asymmetric-top testbed + §6.6 infinite-separation decoupling testbed + §6.7 C2 cross-implementation bit-identity gate; forward-looking OQ-6.1..6.6 for v0.2) | §8.1 free top + §8.2 damped top + §8.3 Lagrange top — all three are analytical testbeds with Podolský-derived closed-form solutions, directly usable as oracles for the validation gates. |

The **v0.2 extension** adds the restricted three-body L2 layer: rotating Earth-Sun frame, pseudo-potential including centrifugal and Coriolis terms, linearised stability analysis at L2, equations of motion for two coupled rigid bodies orbiting in that frame. Reference anchor: Szebehely 1967, Murray & Dermott 1999. Out of scope for v0.1.

## §2. Canonical-tier promise

This documentation tier holds itself to the **lege artis** standard. The same discipline that governs `lege-artis/fourier` (see `fourier/docs/canonical/en/00-introduction.md` §2) applies here:

- Every equation of motion is **derived from first principles** in this tier. We do not say "the gravity-gradient torque is $3 G m / r^3 \cdot \hat{r} \times (\mathbf{I} \hat{r})$ because Hughes 2004 says so" — we show the multipole expansion of mutual gravitational potential, identify the tidal term, derive the corresponding torque from $\boldsymbol{\tau} = -\nabla_\theta U$, and verify the limit reduces correctly.
- Every claim cites a primary source. See `shared/reference-bibliography/refs.bib`. For analytical-mechanics formalism, Thorne & Blandford 2017, Goldstein 2002, Arnold 1989, Marsden & Ratiu 1999, Landau–Lifshitz Vol I. For domain-specific spacecraft attitude dynamics, Hughes 2004. For gravitational potentials of extended bodies, Murray & Dermott 1999, Jackson 1999 (multipole formalism transferred from electrodynamics).
- **Wikipedia is not a primary source.** It can be a sanity check, never an authority.
- **Numerical libraries (NumPy, SciPy, etc.) are oracles, not source material.** When we compare against `scipy.integrate.solve_ivp`, the comparison is for cross-validation against an independent implementation; the canonical equations remain the source of truth.
- The mathematical content lives in **two parallel forms**: prose-with-equations in these `docs/canonical/en/NN-*.md` files, and machine-checkable equation files in `shared/canonical-equations/*.{md,tex}` (queued for v0.1 chapters 02-05).
- Engineer-tier docs (`docs/engineer/en/`, forthcoming) translate this content for practitioners but never contradict it. If a future engineer-tier statement disagrees with this canonical tier, the canonical tier wins; please file an issue.

## §3. Why a canonical tier at all?

The engineer-tier first-cut at `experiments/jwst-l2-first-cut/` is built around the Newton–Euler equations expressed directly in inertial Cartesian coordinates:

$$
\dot{\mathbf{p}}_i = \mathbf{F}_i, \qquad \dot{\mathbf{L}}_i = \boldsymbol{\tau}_i \quad \text{for } i = A, B
$$

with the angular velocity propagated in the body frame via the inertia tensor:

$$
\mathbf{I}\, \dot{\boldsymbol{\omega}} + \boldsymbol{\omega} \times (\mathbf{I}\, \boldsymbol{\omega}) = \boldsymbol{\tau}_{\text{body}}.
$$

This is correct and computes conservation laws to machine precision (the first-cut achieves $|dE/E_0| = 6.2 \times 10^{-12}$ and $|dL/L_0| = 2.1 \times 10^{-10}$ over a 600-second integration at $dt = 0.05$ s). It is, however, *not derived* in the engineer-tier — it is *stated*, and validated against analytical limits (torque-free symmetric top, etc.).

A reader at the canonical tier asks: *why* is the right-hand side of these equations what it is? Where does the angular-velocity-cross-inertia term come from? Why is the gravity-gradient torque $3 G m / r^3 \cdot \hat{r} \times (\mathbf{I} \hat{r})$ to first order in the body-size-over-separation ratio, and what is the next-order term? Why are the Newton–Euler equations correct in inertial coordinates but the Euler equations correct in the body frame?

The canonical tier answers these questions, in order, starting from the variational principle. We treat the configuration of two rigid bodies as a point on the configuration manifold $SE(3) \times SE(3)$ (the special Euclidean group describing rigid motions, twice over), write the Lagrangian $L = T - V$ on its tangent bundle, derive the Euler–Lagrange equations, and specialise. The Newton–Euler equations of the engineer-tier emerge as a particular coordinate expression of this; the Euler equations emerge from Euler–Poincaré reduction by the left-invariance of kinetic energy under body rotations.

The benefit of this approach is not pedagogical alone. Three concrete benefits surface:

1. **Conservation laws fall out structurally.** Energy conservation is invariance under time-translation; total linear and angular momentum conservation are invariances under spatial translation and rotation respectively. Noether's theorem makes these visible at the Lagrangian level, before any numerical integrator touches the problem.
2. **Reduced equations are simpler than reduced-from-Cartesian.** The Euler equations in body frame ($\mathbf{I}\,\dot{\boldsymbol{\omega}} + \boldsymbol{\omega} \times \mathbf{I}\,\boldsymbol{\omega} = \boldsymbol{\tau}$) are a structural consequence of left-invariance of the kinetic energy on $SO(3)$, not a coordinate trick. Once derived this way, they generalise immediately to other Lie groups (e.g., $SE(3)$ for a body that is also translating, with both rotational and translational dynamics intertwined).
3. **Symplectic integrators become obvious.** The engineer-tier uses RK4, which is non-symplectic. The canonical tier's variational formulation suggests Lie-group variational integrators (Hairer–Lubich–Wanner 2006) as the structure-preserving upgrade. We do not implement them in v0.1; we make clear what they would do and why.

## §4. The first-principles requirement, in concrete terms

For each chapter in this tier:

- Every equation has either a derivation in the chapter itself or a citation to a primary source where the derivation appears. No "stated without proof" equations.
- Every claim of physical interpretation is anchored: "this term is the centrifugal pseudo-force experienced in the rotating body frame" must be either derived or cited.
- Numerical examples (computed values, e.g. inertia anisotropies of the JWST-like and probe-like bodies) reference back to the engineer-tier code that produces them, with the line range cited.
- Validation against analytical limits is a chapter requirement (Chapter 06): each canonical equation must be tested against at least one configuration with closed-form analytical solution, and the difference reported.

This discipline mirrors `lege-artis/fourier`'s C1-C4 acceptance criteria (see `_config/LEGE-ARTIS-LANGUAGE-DOCTRINE-v0.1.md` §4.4), specialised to the analytical-mechanics domain. The criteria are:

| Criterion | JWST-L2 interpretation |
|---|---|
| **C1** Textbook-verbatim test coverage | Chapter 06 tests each canonical equation against ≥3 analytical-limit configurations (torque-free symmetric top; binary point-mass Kepler; free precession in central field). **Full catalogue:** `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md` — Class A (exact closed-form solutions) and Class B (edge-case limit configurations), each with mathematical setup, reference anchor, expected numerical tolerance, and test-file location. |
| **C2** Cross-implementation bit-identity (where applicable) | When v0.2 introduces an alternative implementation (e.g., symplectic Lie-group integrator alongside the current RK4), cross-implementation comparison at ULP × √N tolerance becomes a CI requirement. |
| **C3** Per-routine textbook citation in source headers | Engineer-tier code (`dynamics.py`, `geometries.py`) carries header citations to canonical-tier chapter sections, grep-audited pre-build. (Forthcoming in v0.2 of engineer-tier.) |
| **C4** Empirical $O(N \log N)$ complexity audit | Not directly applicable to ODE integration (which is $O(N \cdot N_{steps})$). Substituted: empirical conservation-error vs. step-size scaling audit (expected $O(dt^4)$ for RK4 energy drift). See `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md` §3.6 for the test setup. |

The canonical-tier discipline distinguishes two complementary test classes, both required for honest model validation:

- **Class A — Exact closed-form analytical solutions** as oracles. At any time $t$, the analytical value is ground truth; the numerical value's distance from it is a measurable error. Answers *"is my integrator solving the equations correctly?"*
- **Class B — Edge-case limit configurations** stressing the model's domain of validity. Answers *"am I solving the right equations?"*

Without Class B, a model can be numerically accurate yet physically wrong. Without Class A, the converse. The full catalogue of v0.1 testbeds and edge-cases lives in `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md`; each entry expands into Chapter 06 derivation prose when that chapter lands.

## §5. Reading path

For a reader with strong physics background and prior exposure to Lagrangian mechanics:

1. This Introduction (you are here)
2. Chapter 01 — Configuration manifold and variational principle (next file)
3. Chapter 02 — Kinetic energy on $SE(3) \times SE(3)$
4. Chapter 03 — Mutual gravitational potential and multipole expansion
5. Chapter 04 — Euler–Lagrange equations and Euler–Poincaré reduction
6. Chapter 05 — Body-frame Euler equations as a specialisation
7. Chapter 06 — Analytical limits and validation gates

For a reader without prior Lagrangian-mechanics exposure, start with Feynman Vol II Ch 19 (Principle of Least Action) and Goldstein Ch 1–2, then return here. Thorne & Blandford Vol I provides the same starting material at a higher density.

For a reader interested in the code rather than the math, see the engineer-tier `docs/engineer/en/00-quick-start.md` (forthcoming) or the existing first-cut entry points `experiments/jwst-l2-first-cut/dynamics.py` and `run_first_example.py`. The Shad-tier `SHAD-NARRATIVE.md` is the narrative companion to the first-cut.

## §6. What this tier explicitly does not cover

For honesty and to bound expectations:

- **General relativistic corrections.** JWST orbital velocity $v \approx 0.3 \text{ km/s}$ relative to Earth-Sun barycentre; $v/c \approx 10^{-6}$; rotational coupling effect at first relativistic order scales as $(v/c)^2 \approx 10^{-12}$, which is well below any current measurement precision for rotational state. We work in pure Newtonian gravity.
- **Solar radiation pressure on JWST's sunshield.** This is the dominant non-gravitational perturbation on the actual JWST mission and is responsible for a substantial fraction of the spacecraft's station-keeping fuel budget. It is a real perturbation; we omit it from this canonical tier because (a) it is not a gravitational coupling, (b) it depends on detailed surface optics not available in the canonical equations, and (c) the first-cut engineer-tier prototype focuses on the gravitational coupling specifically. v0.2 of this canonical tier will catalogue what would change if the perturbation were included.
- **Internal vibration modes of either body.** We assume rigid bodies. The actual JWST has structural flexibility; modelling it requires a continuum-mechanics layer outside the rigid-body framework. Out of scope for the lege-artis project.
- **Magnetic and electrostatic interactions.** JWST is uncharged to high precision and operates in a near-vacuum; these are negligible at the L2 environment and ignored without further comment.
- **Long-period perturbations from Moon, Venus, Jupiter.** Mission-design relevance is real; canonical-tier scope confines to the immediate two-body gravitational coupling. Future v0.3 extension would add them.

## §7. Cross-references and companions

- **Companion engineer-tier:** `docs/engineer/en/` (forthcoming v0.2 of engineer docs).
- **Companion Shad-tier:** `SHAD-NARRATIVE.md` in the project root (the engineer-narrated birthday-edition shipment writeup).
- **Reference bibliography:** `shared/reference-bibliography/refs.bib`.
- **Equation files** (forthcoming): `shared/canonical-equations/lagrangian.md`, `…/euler-equations.md`, `…/gravity-gradient-torque.md`, `…/multipole-expansion.md`. Each equation file will be cross-referenced from the chapter section that derives it.
- **First-cut code under derivation here:** `geometries.py`, `dynamics.py`, `run_first_example.py`, `testcases.py`, `tests/`.

---

**Next chapter:** [01 — Configuration manifold and variational principle](01-configuration-manifold.md)

**End 00-introduction.md**
