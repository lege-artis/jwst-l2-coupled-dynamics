# Canonical Validation — Testbeds and Edge-Cases (v0.1)

| Field | Value |
|---|---|
| Document ID | CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1 |
| Status | DRAFT v0.1 — locks validation discipline for canonical-tier Chapter 06; ready for engineer-tier test-suite expansion to cite |
| Authored | 2026-05-23 |
| Project-owner direction | "analytical solutions as testbeds for numerics are important as well as discussion of limit situations to formulate edge-cases for dynamics to be sure output is repeatedly testable and traceable" |
| Sanitisation context | Private monorepo only until first-cut subtree-publish |
| Related | `docs/canonical/en/00-introduction.md` §4 (C1-C4 interpretation); `_config/LEGE-ARTIS-LANGUAGE-DOCTRINE-v0.1.md` §4.4 |

---

## §1 Validation discipline (the contract this spec locks)

For the JWST-L2 canonical-tier work to honestly earn its place in LegeArtis, every numerical integrator output must be **falsifiable** against an externally-anchored reference. Two complementary classes of reference are required, both formally catalogued here:

**Class A — Exact closed-form analytical solutions.** Configurations of the two-rigid-body gravitationally-coupled system for which the equations of motion admit a closed-form integral (constant of motion or explicit time-dependence formula). These serve as **oracles**: at any time $t$, the analytical value is the ground truth, and the numerical value's distance from it is a measurable error.

**Class B — Edge-case limit configurations.** Configurations where the dynamics either (a) reduce to a simpler well-understood problem in some limit (asymptotic regime), or (b) approach a singular boundary of the model's domain of validity. Class B testbeds verify that the integrator handles transitions gracefully and that the model's stated assumptions are respected.

Both classes are required because they answer different questions:
- Class A answers *"is my integrator solving the equations correctly?"* — pure numerical accuracy.
- Class B answers *"am I solving the right equations?"* — model-domain validation.

Each entry below ships with: mathematical setup, what it tests, reference anchor in the bibliography, expected numerical tolerance, and a pointer to the existing or planned test-file location. When Chapter 06 of the canonical tier is fully authored, each entry expands into a derivation block; the catalogue here is the project-stable contract that survives any one chapter revision.

The catalogue covers **v0.1 scope** (two rigid bodies + mutual gravity, no L2 frame). v0.2 will add Class A/B entries for the restricted three-body L2 layer.

---

## §2 Class A — Exact closed-form analytical solutions

### §2.1 — Torque-free symmetric top (single rigid body)

**Mathematical setup.** Single rigid body with axisymmetric inertia tensor $\mathbf{I} = \mathrm{diag}(I_\perp, I_\perp, I_\parallel)$; no external torque ($V \equiv 0$). Euler equations reduce to:

$$
I_\perp \dot{\omega}_x = (I_\perp - I_\parallel)\,\omega_y \omega_z, \quad
I_\perp \dot{\omega}_y = (I_\parallel - I_\perp)\,\omega_x \omega_z, \quad
I_\parallel \dot{\omega}_z = 0.
$$

**Closed-form solution.** $\omega_z(t) = \omega_z(0) = \mathrm{const}$. The transverse components satisfy harmonic oscillation:

$$
\omega_x(t) = \omega_\perp \cos(\lambda t + \phi_0), \qquad
\omega_y(t) = \omega_\perp \sin(\lambda t + \phi_0),
$$

with **Euler precession angular frequency**:

$$
\lambda = \frac{I_\parallel - I_\perp}{I_\perp}\, \omega_z.
$$

Total kinetic energy $E$ and total angular momentum $|\mathbf{L}|$ are exactly conserved.

**What it tests.** Body-frame Euler-equation integration; principal-axis decomposition; precession-frequency extraction from a long-horizon time series. Connects directly to the fourier × JWST-L2 cross-project example (`fourier/examples/jwst-l2-tumble-spectrum/`), where this same precession frequency is recovered from the integrated NDJSON output via `np.fft.rfft`.

**Reference anchor.** Landau-Lifshitz Vol I §35 (closed-form integration of Euler's equations for axisymmetric top).

**Expected numerical tolerance.** RK4 at $dt \le 0.05$ s over $T \le 600$ s: $|dE/E_0| < 10^{-11}$, $|d\mathbf{L}/\mathbf{L}_0| < 10^{-9}$, precession-frequency recovery error $< 0.2$ FFT bins (limited by sampling resolution, not integration accuracy).

**Test-file location.** `tests/test_dynamics.py::TestRK4OnTorqueFreeSymmetricTop` (**implemented in first-cut**; conservation gate passing at $6.2 \times 10^{-12}$ for energy and $2.1 \times 10^{-10}$ for angular momentum).

---

### §2.2 — Torque-free asymmetric top (Dzhanibekov / tennis-racket theorem)

**Mathematical setup.** Single rigid body with three distinct principal moments $I_{xx} < I_{yy} < I_{zz}$; no external torque.

**Closed-form solution.** Rotation about either extreme axis ($I_{xx}$ or $I_{zz}$) is **stable**: small perturbations bounded. Rotation about the middle axis ($I_{yy}$) is **unstable**: small perturbations grow exponentially, the body undergoes a characteristic flip, returns approximately to its initial state, and flips again. The complete solution involves Jacobi elliptic functions; for small perturbations, the instability growth rate is:

$$
\sigma = \omega_y \sqrt{\frac{(I_{yy} - I_{xx})(I_{zz} - I_{yy})}{I_{xx} I_{zz}}}.
$$

**What it tests.** Long-horizon integrator stability across a sequence of flips; correct sign of the instability; conservation of $E$ and $|\mathbf{L}|$ during the chaotic-looking transitions.

**Reference anchor.** Landau-Lifshitz Vol I §37; for the elliptic-function explicit solution, Goldstein §5.6.

**Expected numerical tolerance.** Conservation as in §2.1 over the full integration window including flips. Period between flips reproducible to $< 5\%$ across runs at varying initial perturbation amplitude (the period is amplitude-dependent for the full nonlinear solution).

**Test-file location.** **Not yet implemented.** Queued for v0.2 of engineer-tier test suite; planned as `tests/test_dynamics.py::TestDzhanibekovInstability`.

---

### §2.3 — Two point-mass Kepler problem (monopole-monopole limit)

**Mathematical setup.** Two bodies in the limit where each body's spatial extent is negligible compared to their separation: $|R_A - R_B| \gg \mathrm{body\ sizes}$. Mutual gravitational potential reduces to the monopole-monopole term:

$$
V = -\frac{G m_A m_B}{|\mathbf{R}_A - \mathbf{R}_B|}.
$$

**Closed-form solution.** In the centre-of-mass frame, the relative coordinate $\boldsymbol{\rho} = \mathbf{R}_B - \mathbf{R}_A$ traces a conic section (ellipse, parabola, hyperbola depending on energy). Orbital period for bound (elliptical) orbits:

$$
T_{\mathrm{orbit}} = 2\pi \sqrt{\frac{a^3}{G(m_A + m_B)}}
$$

where $a$ is the semi-major axis. Energy, angular momentum, and Laplace-Runge-Lenz vector are exactly conserved.

**What it tests.** Translational equations of motion in the engineer-tier code; multipole truncation cleanliness (when inertia tensors are set to zero, the rotational coupling vanishes and the system MUST reduce to Kepler).

**Reference anchor.** Landau-Lifshitz Vol I §13-§15.

**Expected numerical tolerance.** Orbital period recovery within $10^{-6}$ relative error over one orbital period; eccentricity stable to similar precision; LRL vector magnitude conserved to $10^{-8}$.

**Test-file location.** **Not yet implemented.** Queued; planned as `tests/test_dynamics.py::TestKeplerLimit` with inertia tensors $\mathbf{I}_A = \mathbf{I}_B = \mathbf{0}$ as the trigger.

---

### §2.4 — Free precession of axisymmetric body in central gravitational field

**Mathematical setup.** Axisymmetric rigid body (body $A$) with $I_{xx,A} = I_{yy,A} \neq I_{zz,A}$, orbiting a spherically symmetric central point mass ($\mathbf{I}_B = \mathbf{0}$, $m_B$ finite). Mutual potential is monopole-monopole only (the point-mass body has no extended-body multipole contribution); gravity-gradient torque on body $A$ vanishes because the potential field due to a point mass is spherically symmetric and exerts no torque on any body whose principal axes are not coupled to position.

**Closed-form solution.** Translational motion: Kelperian orbit as in §2.3. Rotational motion: free Euler precession of body $A$ as in §2.1. The two are **decoupled** — no exchange of angular momentum between orbit and spin.

**What it tests.** That the engineer-tier gravity-gradient-torque implementation correctly evaluates to zero in this limit. Catches sign errors and accidental tidal-from-monopole bugs that would not surface in pure §2.1 or pure §2.3 tests.

**Reference anchor.** Combination of Landau-Lifshitz §13-§15 (orbital) + §35 (rotational); explicit decoupling argument in Goldstein §5.6 / Hughes 2004 Ch 8.

**Expected numerical tolerance.** Cross-coupling between orbital and rotational angular momentum bounded by integration round-off, $< 10^{-10}$ relative over one orbital period.

**Test-file location.** **Not yet implemented.** Planned as `tests/test_integration.py::TestFreeRotationKeplerOrbitDecoupling`.

---

### §2.5 — Small-angle libration about gravity-gradient-stabilised equilibrium

**Mathematical setup.** Axisymmetric body in circular orbit around a central point mass, with body principal axis $\mathbf{e}_z$ in the orbit plane and small perturbations $\theta(t), \phi(t)$ about the radial-pointing equilibrium. To linear order in $(\theta, \phi)$, equations of motion reduce to coupled harmonic oscillators with frequencies determined by the orbit angular velocity $n$ and the inertia ratio.

**Closed-form solution.** Two normal-mode frequencies $\omega_1, \omega_2$ in terms of orbit rate $n$ and inertia ratio; small-amplitude oscillations are sinusoidal.

**What it tests.** Gravity-gradient torque at leading order in body-extent (the multipole-monopole tidal term); small-amplitude regime where linear-ODE theory is reliable.

**Reference anchor.** Hughes 2004 Ch 8 (spacecraft attitude dynamics, gravity-gradient stabilisation).

**Expected numerical tolerance.** Period match to analytical formula within $10^{-4}$ relative error over 10 orbital periods; mode amplitudes stable to $10^{-6}$.

**Test-file location.** **Not yet implemented.** Planned as `tests/test_integration.py::TestGravityGradientLibration`.

---

### §2.6 — Harmonic-oscillator approximation of relative motion near Lagrangian L2 (v0.2 testbed, listed here for completeness)

**Status.** **v0.2 of canonical tier.** Not in v0.1 scope (requires the L2 / restricted three-body extension). Listed here as forward-looking entry; testbed implementation defers until v0.2.

---

## §3 Class B — Edge-case limit configurations

### §3.1 — Probe-as-point-mass limit ($m_B \to 0$, $\mathbf{I}_B \to \mathbf{0}$)

**Setup.** Body B's mass and inertia tensor approach zero.
**Expected behaviour.** System decouples into (a) free torque-free Euler precession of body A in inertial space, and (b) test-particle motion of B in A's gravitational field. Mutual coupling reduces to a one-way force on B by A; no back-reaction.
**What it tests.** Code paths that would silently misbehave in the limit (e.g., $1/m_B$ terms in any equation). Also: the implementation's response to ill-conditioned inertia tensors.
**Tolerance.** Body A's Euler precession matches §2.1; body B's orbit matches §2.3 with $A$ treated as a point mass.
**Test-file location.** **Not yet implemented.** Planned as `tests/test_dynamics.py::TestPointMassLimitForBodyB`.

### §3.2 — Infinite-separation limit ($|\mathbf{R}_A - \mathbf{R}_B| \to \infty$)

**Setup.** Bodies placed at very large separation; mutual potential $V \to 0$.
**Expected behaviour.** Two independent torque-free Euler tops; total energy is exactly the sum of individual rotational kinetic energies; each body's angular momentum is conserved independently.
**What it tests.** Multipole truncation honesty (higher-order multipole terms should decay correctly with $1/r^n$); absence of spurious coupling at large $r$.
**Tolerance.** Cross-coupling between body A's and body B's angular momenta below $10^{-12}$ relative per orbital period.
**Test-file location.** **Not yet implemented.** Planned as `tests/test_dynamics.py::TestInfiniteSeparationDecoupling`.

### §3.3 — Close-approach limit ($|\mathbf{R}_A - \mathbf{R}_B| \to \mathrm{body\ size}$)

**Setup.** Bodies approach distances comparable to their physical sizes.
**Expected behaviour.** **Multipole expansion breaks down.** The current first-cut truncates at quadrupole; physical accuracy degrades smoothly as separation approaches body size. **The code must NOT silently produce wrong answers** — it should either flag the assumption violation or progressively diverge from physical reality in a documented way.
**What it tests.** The implementation's awareness of its own domain-of-validity boundary. Forces explicit documentation of where the model stops being applicable.
**Tolerance.** Quadrupole-truncation error grows as $\sim (\mathrm{body\ size}/r)^4$; when this ratio exceeds $0.3$, an assertion or warning should fire (engineer-tier code, to be added in v0.2).
**Test-file location.** **Not yet implemented.** Planned as `tests/test_dynamics.py::TestCloseApproachAssumptionViolation`.

### §3.4 — Zero-rotation limit ($\boldsymbol{\omega}_A = \boldsymbol{\omega}_B = \mathbf{0}$ at $t=0$)

**Setup.** Both bodies non-rotating at initial time; mutual gravity active.
**Expected behaviour.** Translational dynamics: Kepler orbit (or its degenerate limit if angular momentum is zero) as in §2.3. Rotational dynamics: gravity-gradient torque activates if either body has non-spherical inertia and the principal axes aren't aligned with the local radial direction; otherwise rotation remains zero.
**What it tests.** Initial conditions on a measure-zero subset of phase space; numerical noise should not spuriously activate rotation.
**Tolerance.** When inertia tensors are spherical ($I_{xx} = I_{yy} = I_{zz}$), rotation remains numerically zero ($< 10^{-12}$) through the integration. When inertia tensors are non-spherical, rotation activates from zero and matches the linearised libration prediction §2.5 for small amplitude.
**Test-file location.** **Not yet implemented.** Planned as `tests/test_dynamics.py::TestZeroRotationInitialCondition`.

### §3.5 — Symmetric / aligned initial conditions

**Setup.** Initial $\boldsymbol{\omega}$ exactly aligned with a principal axis of the inertia tensor; or $\Delta\mathbf{R}$ exactly along a principal axis; or both.
**Expected behaviour.** Symmetry is preserved by the equations of motion; the integrator must not introduce numerical drift that breaks the symmetry.
**What it tests.** Time-stepping symmetry preservation; floating-point bit-identity discipline (especially for RK4, which is not symmetric in time and can drift if the integration scheme has hidden asymmetries).
**Tolerance.** Symmetry-violating components stay below $10^{-13}$ over the integration window.
**Test-file location.** **Not yet implemented.** Planned as `tests/test_dynamics.py::TestPrincipalAxisRotationStability` (extension of existing §2.1 test).

### §3.6 — Integrator convergence ($dt \to 0$ and $dt \to \mathrm{large}$)

**Setup.** Run §2.1 setup (torque-free symmetric top) at $dt \in \{0.5, 0.1, 0.05, 0.01, 0.005, 0.001\}$ s.
**Expected behaviour.** Energy-conservation error scales as $O(dt^4)$ for RK4 over a fixed time horizon. At small enough $dt$, error saturates at the floating-point round-off floor ($\sim 10^{-15}$ per step times $N_{\mathrm{steps}}$). At large $dt$, error grows nonlinearly; an instability threshold exists where the integration diverges.
**What it tests.** The numerical method's order-of-accuracy claim (RK4 is 4th-order; this is the empirical confirmation). Also: identifies the practical $dt$ ceiling for the engineer-tier default integration parameters.
**Tolerance.** Log-log fit of $|dE|$ vs $dt$ should have slope $4.0 \pm 0.1$ in the convergent regime.
**Test-file location.** **Not yet implemented.** Planned as `tests/test_dynamics.py::TestRK4OrderOfAccuracy`.

### §3.7 — Long-horizon drift assessment

**Setup.** §2.1 setup integrated over $T = 10^4$ orbital periods at default $dt$.
**Expected behaviour.** RK4 produces a secular drift in energy proportional to $dt^4 \cdot T$ (non-symplectic). The drift is small but cumulative; documented and bounded.
**What it tests.** Whether the engineer-tier integrator is suitable for long-mission-duration simulations (currently it is not — long-horizon should use a symplectic Lie-group integrator per Hairer-Lubich-Wanner 2006).
**Tolerance.** Cumulative $|dE/E_0|$ at $T = 10^4 T_{\mathrm{orbit}}$ should stay below $10^{-4}$ for the default first-cut configuration.
**Test-file location.** **Not yet implemented.** Planned as `tests/test_integration.py::TestLongHorizonRK4Drift`.

---

## §4 Traceability discipline

Every Class A or Class B test entry above ships with the same five fields:
1. **Mathematical setup** — minimum reproducible configuration.
2. **Closed-form solution** (Class A) or **expected behaviour** (Class B) — what the numerical output should match.
3. **Reference anchor** — primary source in `shared/reference-bibliography/refs.bib`.
4. **Expected numerical tolerance** — falsifiable quantitative bound, not "looks reasonable."
5. **Test-file location** — exact pytest function name in the engineer-tier code.

When Chapter 06 of the canonical tier is fully authored, each entry expands into a derivation block: starting from the canonical equations of motion derived in Chapters 02-05, recover the Class A closed-form solution by analytical reduction, or characterise the Class B limit by direct substitution. The test-file then verifies the numerical integrator matches.

This discipline implements language-doctrine §4.4 criterion **C1** (textbook-verbatim test coverage ≥3 cases per published anchor) for the analytical-mechanics domain. **C2** (cross-implementation bit-identity) becomes meaningful only when a symplectic Lie-group integrator is added alongside RK4. **C3** (per-routine textbook citation in source headers) is implemented inside the engineer-tier code's docstring blocks. **C4** (empirical complexity audit) is replaced here by **§3.6 integrator order-of-accuracy** (the analogous "empirical scaling-law audit" for ODE integration).

## §5 v0.1 implementation status — at session close 2026-05-23

| Test ID | Class | Status |
|---|---|---|
| §2.1 Torque-free symmetric top | A | **Implemented** in first-cut |
| §2.2 Dzhanibekov asymmetric top | A | Planned |
| §2.3 Kepler point-mass limit | A | Planned |
| §2.4 Free rotation + Keplerian orbit decoupling | A | Planned |
| §2.5 Gravity-gradient libration | A | Planned |
| §3.1 Point-mass limit for B | B | Planned |
| §3.2 Infinite-separation decoupling | B | Planned |
| §3.3 Close-approach assumption violation | B | Planned (with assertion design) |
| §3.4 Zero-rotation initial conditions | B | Planned |
| §3.5 Principal-axis symmetry preservation | B | Planned (extension of §2.1) |
| §3.6 RK4 order-of-accuracy convergence | B | Planned |
| §3.7 Long-horizon RK4 drift | B | Planned |

**Implementation cost estimate.** Each Class A test: ~30-60 min Python authoring. Each Class B test: ~30-45 min. Total: ~6-9 hours engineer-tier work to bring the test-suite to full coverage of this catalogue. Suitable for delegation to a single comprehensive Sonnet brief when ready, or for incremental authoring across canonical-tier Chapter 06 sessions.

---

## §6 What this spec does not cover (deferred to v0.2)

- **L2-frame restricted three-body testbeds** (linearised stability eigenvalues at L2; halo orbit closed-form approximation; Lyapunov orbits). All deferred to v0.2 of this catalogue alongside the canonical-tier extension.
- **Relativistic 1PN corrections** as a validation testbed. Justified-omitted at v0.1 (per `00-introduction.md` §6); could become a v0.3+ Class A entry if the rigid-body work ever extends into the regime where it matters.
- **Symplectic-integrator implementation and cross-implementation bit-identity tests** (doctrine C2 in full). Activates only when a second integrator is implemented alongside RK4.
- **Engineer-tier code that fires `assert` statements on §3.3 close-approach violations.** Currently the first-cut silently progresses with degraded physics; the assertion design is queued.

---

**Cross-references**
- `docs/canonical/en/00-introduction.md` §4 — canonical-tier interpretation of C1-C4
- `docs/canonical/en/01-configuration-manifold.md` — geometric foundation for the §6 chapter
- `shared/reference-bibliography/refs.bib` — primary-source anchors
- `_config/LEGE-ARTIS-LANGUAGE-DOCTRINE-v0.1.md` §4.4 — workspace-wide acceptance criteria
- Future `docs/canonical/en/06-validation.md` — chapter where each entry above expands into derivation prose

---

**End CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md**
