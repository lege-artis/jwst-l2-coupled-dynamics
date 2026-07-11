# Chapter 6 — Analytical limits and validation gates

> **Prerequisite reading.** Chapters 01–05. This chapter operates on the boxed total equation set of Chapter 04 §4.7 + the principal-axis specialisation of Chapter 05, validating each piece against analytical-limit configurations catalogued in `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md`.
>
> **Status.** v0.1 (canonical-tier round 3, 2026-05-25). Scaffold + content for §6.1–§6.7 testbeds; §6.8 establishment + §6.9 OQ list for round 4+ work.

---

## §6.1 What this chapter does

The previous five chapters constructed the canonical-tier equations of motion from first principles: configuration manifold + variational principle (Ch 01) → kinetic energy as left-invariant Riemannian metric (Ch 02) → mutual gravitational potential via multipole expansion (Ch 03) → Lagrangian assembly + Euler-Poincaré reduction (Ch 04) → principal-axis Euler equations (Ch 05). The result is the boxed Ch 04 §4.7 equation set on the 18-dimensional reduced phase space $T^* Q_{\mathrm{CoM}}$, with rotational dynamics specialised to the principal-axis form in Ch 05 §5.4 boxed.

A canonical-tier derivation that is internally consistent but unfalsifiable would not earn its place under the lege-artis discipline. Chapter 06 brings the derived equations into contact with **falsifiable analytical-limit testbeds**: configurations of the system for which the equations of motion admit closed-form analytical solutions, against which any numerical implementation can be quantitatively compared. The discipline is the doctrine §4.4 C1 criterion (textbook-verbatim test coverage ≥3 cases per published anchor), adapted to the analytical-mechanics domain.

Following the convention locked in `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md` §1, the testbeds split into two classes:

- **Class A — Exact closed-form analytical solutions.** Configurations where the equations of motion reduce to a system whose solution is known in closed form. The analytical solution serves as an **oracle**: at any time $t$, the analytical value is ground truth, the numerical value's distance from it is a measurable error. Class A answers *"is the integrator solving the equations correctly?"*
- **Class B — Edge-case limit configurations.** Configurations where the dynamics either reduce to a simpler problem in some asymptotic regime, or approach a singular boundary of the model's domain of validity. Class B testbeds verify that the integrator handles transitions gracefully and that the model's stated assumptions are respected. Class B answers *"am I solving the right equations?"*

Both classes are required because they answer different questions. Without Class B, a model can be numerically accurate yet physically wrong. Without Class A, the converse.

Each testbed below ships with: mathematical setup, derivation of the analytical-limit reduction or expected behaviour, closed-form solution (Class A) or expected behaviour (Class B), Podolský 2018 CS-language anchor where available, expected numerical tolerance, and pointer to the engineer-tier / Phase 1 Fortran test-file location. The v0.1 catalogue catalogue is structured around the **six core testbeds** §6.2–§6.6 plus the procedural §6.7 cross-implementation gate. v0.2 will extend with the L2 restricted-three-body family.

## §6.2 Kepler-limit testbed (Class A)

### §6.2.1 Mathematical setup

Set both inertia tensors to zero: $\mathbf{I}_A \to \mathbf{0}, \mathbf{I}_B \to \mathbf{0}$. By Ch 03 §3.4.3 boxed, both tidal pieces $V_{\mathrm{tidal},A}, V_{\mathrm{tidal},B}$ vanish identically (each is proportional to the inertia trace + bilinear inertia projection). The full mutual potential collapses to the monopole-monopole Kepler term:

$$
V \;\xrightarrow{\mathbf{I}_{A,B} \to \mathbf{0}}\; -\frac{G m_A m_B}{|\boldsymbol{\rho}|}.
$$

By Ch 04 §4.7.2 numerical sanity check, the rotational equations decouple (both $\boldsymbol{\tau}_X^{\mathrm{body}} = \mathbf{0}$ and the back-reaction force $\mathbf{F}_{\mathrm{back-reaction},\, X} = \mathbf{0}$ vanish because both scale with $\mathbf{I}_X$). The boxed §4.7 equation set reduces to:

$$
\mu \ddot{\boldsymbol{\rho}} \;=\; -\frac{Gm_A m_B}{|\boldsymbol{\rho}|^3}\, \boldsymbol{\rho}, \qquad \dot{\boldsymbol{\Omega}}_X = \mathbf{0}.
$$

The translational equation is the standard reduced-mass Kepler problem; rotational equations decouple to $\boldsymbol{\Omega}_X = \mathrm{const}$.

### §6.2.2 Analytical solution

In the centre-of-mass frame, the relative coordinate $\boldsymbol{\rho}$ traces a conic section (ellipse, parabola, hyperbola) determined by the energy $E$ and angular momentum $|\mathbf{L}|$ of the relative motion. For a bound (elliptical) orbit with semi-major axis $a$, the orbital period is:

$$
T_{\mathrm{orbit}} \;=\; 2\pi \sqrt{\frac{a^3}{G(m_A + m_B)}},
$$

(Kepler's third law). The orbital energy is $E = -G m_A m_B / (2a)$, the eccentricity is $e = \sqrt{1 + 2 E |\mathbf{L}|^2 / (\mu (Gm_Am_B)^2)}$, and the Laplace-Runge-Lenz vector $\mathbf{A} = \mathbf{p}_{\mathrm{rel}} \times \mathbf{L} - \mu Gm_Am_B \hat{\boldsymbol{\rho}}$ is exactly conserved (its conservation reflects the hidden $SO(4)$ symmetry of the Kepler problem; see Landau-Lifshitz Vol I §15).

### §6.2.3 Podolský anchor (CS-language)

The Kepler problem is treated in classical Czech-language analytical-mechanics textbooks but is not specifically anchored in Podolský 2018 Kapitola 6–8 (which focus on rigid-body kinematics + dynamics + applications). For CS readers, the standard precedent is Landau-Lifshitz Vol I §13–§15 (the *Mechanika* CS translation is widely available in Czech university libraries).

### §6.2.4 Numerical tolerance + test-file

| Quantity | Tolerance | Test-file |
|---|---|---|
| Orbital period vs $2\pi\sqrt{a^3/(GM)}$ | $< 10^{-6}$ relative | `tests/test_dynamics.py::TestKeplerLimit` (planned) |
| Eccentricity stability over $T_{\mathrm{orbit}}$ | $< 10^{-6}$ | (same) |
| LRL vector magnitude over $T_{\mathrm{orbit}}$ | $< 10^{-8}$ | (same) |
| Conservation of $E$ and $|\mathbf{L}|$ | $< 10^{-10}$ | (same) |

Implementation status: planned (per `_specs/...` §2.3). Engineer-tier path: set `geometries.make_zero_inertia_body()` (synthetic constructor, queued) to produce zero-inertia placeholders; invoke `dynamics.py::integrate` with the otherwise-standard configuration; compare against `scipy.integrate.solve_ivp` Kepler-equation reference at the same $dt$.

## §6.3 Free symmetric top testbed (Class A)

### §6.3.1 Mathematical setup

Single rigid body with axisymmetric inertia tensor $\mathbf{I} = \mathrm{diag}(I_\perp, I_\perp, I_\parallel)$; no external torque ($\boldsymbol{\tau}^{\mathrm{body}} = \mathbf{0}$). The boxed Ch 05 §5.2.1 Euler equations specialise (per Ch 05 §5.3.2) to:

$$
I_\perp \dot\Omega_x = (I_\perp - I_\parallel)\,\Omega_y \Omega_z, \quad
I_\perp \dot\Omega_y = -(I_\perp - I_\parallel)\,\Omega_x \Omega_z, \quad
I_\parallel \dot\Omega_z = 0.
$$

The torque-free case is achievable in our two-body setting by (a) setting the other body's mass to zero so there is no Kepler partner and no gravity-gradient field, or (b) placing the other body at infinite separation so the tidal coupling vanishes (Ch 03 §3.8.4 / §6.6 below). Implementationally, the cleanest engineer-tier path is (a) with a single-body integration entry-point, which is what the Phase 1 D8 testbed `backends/fortran/tests/test_symmetric_top.f90` implements.

### §6.3.2 Analytical solution

$\Omega_z(t) = \Omega_z(0) = \mathrm{const}$ (conserved). The transverse components satisfy coupled harmonic oscillation:

$$
\Omega_x(t) = A \cos(\lambda t + \phi_0), \qquad
\Omega_y(t) = A \sin(\lambda t + \phi_0), \qquad
\Omega_z(t) = \Omega_z^0,
$$

with the **Euler precession angular frequency** (often called the "body-cone half-angle precession rate"):

$$
\lambda \;=\; \frac{I_\parallel - I_\perp}{I_\perp}\, \Omega_z^0.
$$

The body-frame angular-velocity vector $\boldsymbol{\Omega}(t)$ traces a **body cone** of half-angle $\arctan(A/\Omega_z^0)$ around the symmetry axis $\hat{\mathbf{e}}_3$ at rate $\lambda$. Conserved quantities: kinetic energy $T = \tfrac{1}{2}(I_\perp A^2 + I_\parallel \Omega_z^{0\,2})$ and angular-momentum magnitude $|\mathbf{L}|^2 = I_\perp^2 A^2 + I_\parallel^2 \Omega_z^{0\,2}$ both exactly conserved.

### §6.3.3 Podolský anchor (CS-language)

**Direct anchor.** Podolský 2018 §8.1 *Volný symetrický setrvačník* derives this exact solution from his boxed Eq. (7.15) torque-free + axisymmetric inertia. His Eq. (8.1) is the analog of our $(\Omega_x, \Omega_y, \Omega_z)$ solution above, and his Eq. (8.2) gives the body-cone geometric interpretation. The Phase 1 D8 Fortran testbed `backends/fortran/tests/test_symmetric_top.f90` uses Podolský Eqs. (8.1)–(8.2) as the CS-language analytical oracle in the comparison.

### §6.3.4 Numerical tolerance + test-file

| Quantity | Tolerance | Test-file |
|---|---|---|
| $\Omega_z$ drift over $T = 600$ s | $< 10^{-12}$ relative | `backends/fortran/tests/test_symmetric_top.f90` (Phase 1 D8) + `tests/test_dynamics.py::TestRK4OnTorqueFreeSymmetricTop` (engineer-tier, implemented in first-cut) |
| $A = |\Omega_x(t) + i\Omega_y(t)|$ envelope stability | $< 10^{-12}$ relative | (same) |
| Precession frequency $\lambda$ recovery error | $< 0.2$ FFT bins (limited by sampling, not integration) | (same) |
| Energy / angular-momentum conservation | $|dE/E_0| < 10^{-11}$, $|d\mathbf{L}/L_0| < 10^{-9}$ | (same) |

**Status: implemented in first-cut** at engineer-tier (the `TestRK4OnTorqueFreeSymmetricTop` conservation gate currently passing at $6.2 \times 10^{-12}$ for energy and $2.1 \times 10^{-10}$ for angular momentum; Ch 00 §3 cites these as the headline first-cut diagnostics). Phase 1 Fortran cross-port D8 is the lege-artis bit-identity cross-check.

## §6.4 Gravity-gradient libration testbed (Class A linearised + Class A nonlinear via Lagrange-top reduction)

### §6.4.1 Mathematical setup

Axisymmetric body $A$ ($I_{xx,A} = I_{yy,A} \equiv I_\perp$, $I_{zz,A} \equiv I_\parallel$) in a circular orbit around a point-mass body $B$ ($\mathbf{I}_B = \mathbf{0}$), with the body's symmetry axis $\hat{\mathbf{e}}_3^A$ nominally aligned with $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$ (radial-pointing equilibrium). The boxed Ch 05 §5.4 substituted Euler equations specialise with $\boldsymbol{\tau}_A^{\mathrm{body}}$ from gravity-gradient + axisymmetric $\mathbf{I}_A$. The orbit rate is $n = \sqrt{G(m_A + m_B)/|\boldsymbol{\rho}|^3}$.

Perturb the orientation by small angles $(\theta, \phi)$ about the equilibrium (a pitch–roll perturbation, with $\theta$ rotation about body $\hat{\mathbf{e}}_1$ and $\phi$ about body $\hat{\mathbf{e}}_2$). To linear order in $(\theta, \phi)$, the body-frame separation direction becomes $\hat{\boldsymbol{\rho}}_{\mathrm{body},A} \approx (-\phi, -\theta, 1)$, and the Ch 05 §5.4 substituted equations give the linearised libration equations (Ch 05 §5.4.1 derivation extended to both transverse modes; Ch 03 §3.7.2 boxed orientation Jacobian provides the equivalent matrix form).

### §6.4.2 Closed-form solution — linearised regime

The linearised libration is two-dimensional, with normal-mode frequencies determined by the inertia ratio $(I_\parallel - I_\perp)/I_\perp$ and the orbit rate $n$. For the planar pitch mode about the radial-pointing equilibrium (Ch 05 §5.4.1 derivation; cross-checked against Ch 03 §3.7.2 orientation Jacobian):

$$
\omega_{\mathrm{lib}}^2 \;=\; \frac{3 G m_B (I_\perp - I_\parallel)}{I_\perp\, |\boldsymbol{\rho}|^3}
\;=\; 3 n^2\, \frac{I_\perp - I_\parallel}{I_\perp},
$$

using $n^2 = G m_B / |\boldsymbol{\rho}|^3$ for the orbit rate in the point-mass-primary limit (Ch 06 §6.2 Kepler limit). For prolate ($I_\parallel < I_\perp$): $\omega_{\mathrm{lib}}^2 > 0$ → stable libration. For oblate ($I_\parallel > I_\perp$): $\omega_{\mathrm{lib}}^2 < 0$ → exponential instability; the symmetry-axis-along-radial attitude is the unstable equilibrium, and the body relaxes toward the symmetry-axis-perpendicular-to-radial attitude. The full pitch-roll-yaw three-mode coupling treatment is in Hughes 2004 §8.4 Eqs. (8.40)–(8.46).

### §6.4.3 Closed-form solution — full nonlinear via Lagrange-top effective potential (resolves OQ-5.6 + OQ-6.4)

The finite-amplitude libration about the radial-pointing equilibrium is the spacecraft-attitude analog of the classical **Lagrange top** problem (heavy symmetric top with fixed point in uniform gravity). The connection is structural: an axisymmetric body in a position-dependent gravitational field with the symmetry axis tilted at angle $\vartheta$ from the field direction satisfies a Lagrangian with the same algebraic form regardless of whether the "field" is a uniform gravitational acceleration (Podolský §8.3) or a position-dependent gravity-gradient torque (our §6.4 setup). The Podolský 2018 §8.3 derivation transfers verbatim with substitution $g \to 3 G m_B / |\boldsymbol{\rho}|^3$ as the effective-acceleration parameter, and the Lagrange-top **effective-potential method** provides the finite-amplitude closed-form solution.

**Setup.** Adopt the intrinsic-ZXZ Euler-angle chart of Ch 04 §4.3.1 / Podolský §6.4 Eqs. (6.14)–(6.16) with body-$A$ symmetry axis along body $\hat{\mathbf{e}}_3$. At the radial-pointing equilibrium, the Euler angles are $(\varphi, \vartheta, \psi) = (\varphi_0(t), 0, \psi_0(t))$, with $\varphi$ encoding the orbital motion and $\psi$ the spin about the symmetry axis. Perturbation off the equilibrium gives $\vartheta \neq 0$ — the **nutation angle** between body-$\hat{\mathbf{e}}_3$ and $\hat{\boldsymbol{\rho}}$.

The Lagrangian for the axisymmetric body in the gravity-gradient field (Ch 04 §4.2 boxed, specialised to axisymmetric inertia and Euler-angle chart of Ch 04 §4.3) reads:

$$
L \;=\; \tfrac{1}{2} I_\perp (\dot\vartheta^2 + \dot\varphi^2 \sin^2\vartheta)
\;+\; \tfrac{1}{2} I_\parallel (\dot\psi + \dot\varphi \cos\vartheta)^2
\;-\; V_{\mathrm{tidal},\, A}(\vartheta).
$$

In the regime where the orbital radius is large compared to body extent and the orbit period much longer than the spin period — both well-satisfied at hypothetical L2-like configurations (Ch 03 §3.7.2 multi-scale time-scale table) — the tidal potential becomes effectively orientation-dependent only through $\vartheta$:

$$
V_{\mathrm{tidal},\, A}(\vartheta)
\;=\;
-\frac{G m_B}{2 |\boldsymbol{\rho}|^3}\,
\big[\, \mathrm{tr}(\mathbf{I}_A) - 3 (I_\perp \sin^2\vartheta + I_\parallel \cos^2\vartheta)\,\big]
\;=\;
-\frac{G m_B}{2 |\boldsymbol{\rho}|^3}\,
\big[\, (I_\perp + I_\parallel) - 3 I_\perp + 3 (I_\perp - I_\parallel) \cos^2\vartheta \,\big],
$$

up to a $\vartheta$-independent constant absorbed into the energy reference. The orientation-dependent piece is **$\propto \cos^2\vartheta$**, which is precisely the algebraic form of Podolský's gravitational potential in the Lagrange-top setup (his Eqs. 8.7–8.9 with the substitution noted above). The Lagrangian therefore has **cyclic coordinates** $\varphi$ and $\psi$ — neither appears explicitly in $L$ — yielding two first integrals via Noether's theorem (Ch 02 §2.6.1 + Ch 04 §4.6.4 conservation framework):

$$
L_\varphi \;\equiv\; \frac{\partial L}{\partial \dot\varphi}
\;=\; I_\perp \dot\varphi \sin^2\vartheta + I_\parallel (\dot\psi + \dot\varphi \cos\vartheta)\cos\vartheta
\;=\; \mathrm{const},
$$

$$
L_\psi \;\equiv\; \frac{\partial L}{\partial \dot\psi}
\;=\; I_\parallel (\dot\psi + \dot\varphi \cos\vartheta)
\;=\; \mathrm{const}.
$$

These are the **azimuthal angular-momentum component** (about $\hat{\boldsymbol{\rho}}$) and the **body-frame spin angular momentum** (about body-$\hat{\mathbf{e}}_3$), respectively — exactly Podolský's Eqs. (8.10)–(8.11). With $L_\psi$ + $L_\varphi$ conserved, the equation of motion for the remaining degree of freedom $\vartheta$ reduces by quadrature: substitute $\dot\psi$ and $\dot\varphi$ from the first integrals back into the total mechanical energy $E = T + V$ (the third conserved quantity), and solve for $\dot\vartheta^2$ as a function of $\vartheta$ alone.

**Effective potential.** Carry through the substitution. From $L_\psi = I_\parallel(\dot\psi + \dot\varphi\cos\vartheta)$, we have $\dot\psi + \dot\varphi\cos\vartheta = L_\psi/I_\parallel$, and so the spin kinetic energy $\tfrac{1}{2}I_\parallel(\dot\psi + \dot\varphi\cos\vartheta)^2 = L_\psi^2/(2I_\parallel) = \mathrm{const}$ — absorbed into the energy reference. From $L_\varphi = I_\perp\dot\varphi\sin^2\vartheta + L_\psi\cos\vartheta$, we solve $\dot\varphi = (L_\varphi - L_\psi\cos\vartheta)/(I_\perp\sin^2\vartheta)$. Substituting into the rotational kinetic energy:

$$
T_{\mathrm{rot,\, effective}}
\;=\;
\tfrac{1}{2} I_\perp \dot\vartheta^2
\;+\;
\frac{(L_\varphi - L_\psi \cos\vartheta)^2}{2 I_\perp \sin^2\vartheta}
\;+\; \mathrm{const}.
$$

The first term is the kinetic energy of the $\vartheta$ degree of freedom; the second is the **centrifugal-barrier-like** contribution from the eliminated cyclic coordinates. Defining the **effective potential** $V_{\mathrm{ef}}(\vartheta)$ that combines this centrifugal term with the orientation-dependent gravity-gradient piece:

$$
\boxed{
V_{\mathrm{ef}}(\vartheta)
\;=\;
\frac{(L_\varphi - L_\psi \cos\vartheta)^2}{2 I_\perp \sin^2\vartheta}
\;-\;
\frac{3 G m_B (I_\perp - I_\parallel)}{2 |\boldsymbol{\rho}|^3}\, \cos^2\vartheta.
}
$$

The full one-dimensional equation of motion in $\vartheta$ is then:

$$
\tfrac{1}{2} I_\perp \dot\vartheta^2 \;+\; V_{\mathrm{ef}}(\vartheta) \;=\; E_{\mathrm{ef}} \;=\; \mathrm{const}.
$$

This is the **Lagrange-top reduction**: a two-degree-of-freedom rotational problem (after the spin is integrated out) reduced to a single one-dimensional motion of $\vartheta$ in an effective potential. Podolský's Eqs. (8.12)–(8.15) give exactly this structure with our substitution.

**Closed-form solution by quadrature.** The trajectory $\vartheta(t)$ is obtained by separation of variables:

$$
\boxed{
t \;-\; t_0 \;=\;
\int_{\vartheta_0}^{\vartheta(t)}\, \frac{d\vartheta'}{\sqrt{2(E_{\mathrm{ef}} - V_{\mathrm{ef}}(\vartheta'))/I_\perp}}.
}
$$

The integral is **elliptic** in general (the cubic algebraic equation $E_{\mathrm{ef}} = V_{\mathrm{ef}}(\vartheta)$ for the turning points has three roots $\vartheta_1, \vartheta_2, \vartheta_3$; the integral reduces to elliptic integrals of the first and third kinds, expressed in terms of Jacobi elliptic functions or equivalently the Weierstrass $\wp$ function). At small amplitude — $|\vartheta| \ll 1$ — the elliptic integral reduces to the linearised libration of §6.4.2: $\vartheta(t) = \vartheta_{\mathrm{amp}} \cos(\omega_{\mathrm{lib}} t + \phi_0)$ with the linearised frequency formula recovered exactly. At finite amplitude, the period **lengthens** with amplitude (period-amplitude scaling: $T_{\mathrm{lib}}(\vartheta_{\mathrm{amp}})/T_{\mathrm{lib}}(0) = (2/\pi) K(k)$ where $K$ is the complete elliptic integral of the first kind and $k$ is the elliptic modulus determined by $\vartheta_{\mathrm{amp}}$, the spin $L_\psi$, and the inertia ratio).

**Three precession regimes (Podolský §8.3 (a)/(b)/(c)).** Depending on the ratio $L_\psi / L_\varphi$ and the energy $E_{\mathrm{ef}}$, the effective-potential well admits three qualitatively distinct trajectories on the unit sphere (the body-$\hat{\mathbf{e}}_3$ tip):

- **Regime (a) — pure precession:** $\vartheta = \mathrm{const}$ (the body-$\hat{\mathbf{e}}_3$ traces a circle on the sphere). The condition is $dV_{\mathrm{ef}}/d\vartheta = 0$ at the chosen $\vartheta$ value; physically, the spin-induced centrifugal-barrier exactly balances the gravity-gradient torque.
- **Regime (b) — nutation between two latitudes $\vartheta_1 < \vartheta_2$:** the standard libration regime; $\vartheta$ oscillates between turning points $\vartheta_1$ and $\vartheta_2$ where $V_{\mathrm{ef}}(\vartheta_i) = E_{\mathrm{ef}}$. The body-$\hat{\mathbf{e}}_3$ traces a wavy path between two parallels on the unit sphere; cusps appear at the turning points when $\dot\varphi = 0$ there.
- **Regime (c) — looping nutation:** if $L_\psi < L_\varphi/(1+|\cos\vartheta_{\max}|)$ or similar threshold condition, the azimuthal angular velocity $\dot\varphi$ changes sign within the nutation cycle; the trajectory shows loops on the sphere.

The numerical-amplitude scaling of the libration period $T_{\mathrm{lib}}(\vartheta_{\mathrm{amp}})$ — the **finite-amplitude correction** to the linearised §6.4.2 formula — is the key experimentally-falsifiable prediction of the closed-form solution, with explicit dependence on the elliptic modulus.

### §6.4.4 Podolský anchor (CS-language) — Lagrange-top variant

**Direct anchor.** Podolský 2018 §8.3 *Těžký symetrický setrvačník s pevným bodem* (Lagrange top) provides the CS-language analytical oracle. Specifically:

- Podolský Eqs. (8.10)–(8.11) are our cyclic-coordinate first integrals $L_\psi, L_\varphi$ — identical algebraic structure with $g \to 3Gm_B/|\boldsymbol{\rho}|^3$ substitution.
- Podolský boxed Eq. (8.15) is our §6.4.3 boxed effective potential $V_{\mathrm{ef}}(\vartheta)$ — the centrifugal-barrier term identical; the gravitational term substitutes uniform-gravity $V = mgl\cos\vartheta$ with our $-3Gm_B(I_\perp-I_\parallel)\cos^2\vartheta/(2|\boldsymbol{\rho}|^3)$. The functional form differs (his $\cos\vartheta$ vs our $\cos^2\vartheta$) because Podolský's setup has the body's *centre of mass* displaced from the fixed point (giving the cosine dependence on tilt angle from gravity), while ours has a body *tilted* from a radial direction in a position-dependent gravity field (giving cosine-squared dependence). The qualitative phenomenology — bounded nutation, three regimes (a)/(b)/(c) — is identical.
- Podolský Eqs. (8.13) with cases (a)/(b)/(c) maps directly onto our §6.4.3 three regimes.

The Phase 2 testbed `backends/fortran/tests/test_libration.f90` (queued for Phase 2 of the Fortran lege-artis port) will cross-validate against the boxed §6.4.3 $V_{\mathrm{ef}}$ formula numerically: integrate the boxed Ch 04 §4.7 equation set with axisymmetric body in the gravity-gradient field, extract the libration period at multiple amplitudes, and compare against the analytical period $T_{\mathrm{lib}}(\vartheta_{\mathrm{amp}}) = (2/\pi) K(k(\vartheta_{\mathrm{amp}})) \cdot T_{\mathrm{lib}}^{\mathrm{linear}}$. **(Resolves OQ-5.6 + OQ-6.4 jointly.)**

### §6.4.5 Numerical tolerance + test-file

| Quantity | Tolerance | Test-file |
|---|---|---|
| Linearised libration period vs $2\pi/\omega_{\mathrm{lib}}$ | $< 10^{-4}$ relative over 10 orbital periods | `backends/fortran/tests/test_libration.f90` (Phase 2; queued) + `tests/test_integration.py::TestGravityGradientLibration` (planned engineer-tier) |
| Finite-amplitude libration period vs $\boxed{T_{\mathrm{lib}}(\vartheta_{\mathrm{amp}}) = (2/\pi) K(k) \cdot T_{\mathrm{lib}}^{\mathrm{linear}}}$ (§6.4.3) | $< 10^{-3}$ relative at amplitudes up to $\vartheta_{\mathrm{amp}} \lesssim \pi/4$; degrades systematically at larger amplitudes as elliptic-modulus expansion breaks down | (same; Phase 2 finite-amplitude scan) |
| Amplitude stability (no spurious damping or growth) | $< 10^{-6}$ relative over 10 orbital periods | (same) |
| Stability classification (prolate-stable vs oblate-unstable) | Sign-correct from D5 orientation Jacobian eigenvalue (regression guard) | `backends/fortran/tests/test_orientation_jacobian.f90` (Phase 1 D5) |
| Three-regime classification (a)/(b)/(c) | Trajectory on body-$\hat{\mathbf{e}}_3$-tip unit sphere matches Podolský §8.3 phenomenology at corresponding $L_\psi / L_\varphi$ ratios | Phase 2 supplementary |

With §6.4.3's closed-form solution, the full testbed is Class A (both linearised regime and finite-amplitude regime now have analytical oracles). The previous "Class B numerical-only" designation for the finite-amplitude regime is **retired** — the Lagrange-top reduction provides the missing analytical content.

## §6.5 Dzhanibekov asymmetric-top testbed (Class A linearised + Class A nonlinear via Jacobi-elliptic functions)

### §6.5.1 Mathematical setup

Single rigid body with three distinct principal moments $I_{xx} < I_{yy} < I_{zz}$ ; no external torque. Set initial angular velocity nominally about the intermediate axis $\hat{\mathbf{e}}_2$: $\boldsymbol{\Omega}(0) = (0, \Omega_2^0, 0)$ + small perturbation $(\delta\Omega_x, \delta\Omega_y, \delta\Omega_z)$ with $|\delta\boldsymbol{\Omega}| \ll |\Omega_2^0|$. The boxed Ch 05 §5.2.1 Euler equations linearised about this initial state give (Ch 05 §5.3.3 derivation):

$$
I_{xx} \delta\dot\Omega_x = (I_{yy} - I_{zz})\, \delta\Omega_z\, \Omega_2^0, \qquad
I_{zz} \delta\dot\Omega_z = (I_{xx} - I_{yy})\, \delta\Omega_x\, \Omega_2^0.
$$

(The $\delta\Omega_y$ equation is decoupled at linear order; the perturbation initially conserves $\delta\Omega_y \approx 0$ until nonlinear regime.) Substituting one equation into the other:

$$
\delta\ddot\Omega_x \;=\; -\frac{(I_{yy} - I_{xx})(I_{zz} - I_{yy})}{I_{xx} I_{zz}}\, (\Omega_2^0)^2\, \delta\Omega_x.
$$

The coefficient on the right is **negative** because $(I_{yy} - I_{xx}) > 0$ and $(I_{zz} - I_{yy}) > 0$ for the assumed ordering, so the bracketed numerator is positive and the leading sign makes the equation $\delta\ddot\Omega_x = +\sigma^2 \delta\Omega_x$ with $\sigma^2 > 0$. The solution is exponentially growing: $\delta\Omega_x \sim e^{\sigma t}$.

### §6.5.2 Closed-form linearised growth rate

The linearised growth rate is:

$$
\sigma \;=\; \Omega_2^0\, \sqrt{\frac{(I_{yy} - I_{xx})(I_{zz} - I_{yy})}{I_{xx} I_{zz}}}.
$$

This characterises the early-time exponential growth $|\delta\Omega_x(t)| \sim e^{\sigma t}$ before nonlinearities saturate the perturbation. The cross-validation with the textbook tennis-racket-theorem classification (Ch 05 §5.2.2) is exact.

### §6.5.3 Closed-form finite-amplitude solution via Jacobi elliptic functions (resolves OQ-5.5)

For the **finite-amplitude** dynamics — the full Dzhanibekov flip-flip-flip pattern with non-infinitesimal perturbation amplitude — the torque-free asymmetric-top problem is **integrable by quadrature** because the boxed Ch 05 §5.2.1 Euler equations admit *two* conserved quantities (Ch 05 §5.5.1): the body-frame kinetic energy $T$ and the angular-momentum-squared $|\mathbf{L}|^2$. Together with the body-frame Euler equations, these two constraints reduce the three-dimensional system to one effective degree of freedom, whose evolution is given in closed form by **Jacobi elliptic functions**.

**Reduction via conservation laws.** Denote $L_i \equiv I_i \Omega_i$ (the body-frame components of angular momentum; not to be confused with the inertial-frame total angular momentum). The two conserved quantities are:

$$
2T \;=\; L_x^2/I_1 + L_y^2/I_2 + L_z^2/I_3 \;=\; \mathrm{const}, \qquad
|\mathbf{L}|^2 \;=\; L_x^2 + L_y^2 + L_z^2 \;=\; \mathrm{const}.
$$

In body-frame $\mathbf{L}$-space, the trajectory lies on the intersection of:
- The **angular-momentum sphere** $L_x^2 + L_y^2 + L_z^2 = |\mathbf{L}|^2$ (a sphere of radius $|\mathbf{L}|$).
- The **kinetic-energy ellipsoid** $L_x^2/I_1 + L_y^2/I_2 + L_z^2/I_3 = 2T$ (an ellipsoid with semi-axes $\sqrt{2 T I_i}$ along the principal axes).

Generically (for $|\mathbf{L}|^2 \neq 2 T I_i$ for any $i$), the intersection is a smooth closed curve. The trajectory traverses this curve periodically — the **polhode** in body-frame $\mathbf{L}$-space (Goldstein 2002 §5.6 nomenclature).

**Parametrisation by elliptic functions.** Assume $I_1 < I_2 < I_3$ and $|\mathbf{L}|^2 > 2 T I_2$ (the branch whose polhode encircles the **largest**-moment principal axis $\hat{\mathbf{e}}_3$ — this is the branch realised by the §6.5.5 Dzhanibekov test configurations, e.g. $\mathbf{I} = \mathrm{diag}(100,200,400)$, $\boldsymbol{\Omega}_0 = (10^{-3}, 10^{-2}, 10^{-3})$ gives $|\mathbf{L}|^2/(2TI_2) = 1.017$; the complementary branch $|\mathbf{L}|^2 < 2 T I_2$ encircles the smallest-moment axis $\hat{\mathbf{e}}_1$ and follows by swapping the roles of indices 1 and 3 throughout). Substituting the conservation laws into the Euler equation for $L_y$ and using $dL_y/dt = (I_3 - I_1)(L_z L_x)/(I_3 I_1)$, after algebra:

$$
\boxed{
\left(\frac{dL_y}{dt}\right)^2
\;=\;
\frac{(I_2 - I_1)(I_3 - I_2)}{I_1 I_2^2 I_3}\,
\big[\, \beta^2 - L_y^2 \,\big]\,
\big[\, \gamma^2 - L_y^2 \,\big],
}
$$

where $\beta^2 = I_2\,(2 T I_3 - |\mathbf{L}|^2)/(I_3 - I_2)$ and $\gamma^2 = I_2\,(|\mathbf{L}|^2 - 2 T I_1)/(I_2 - I_1)$ — set by the two conserved quantities. On this branch $\beta^2 \le \gamma^2$ (equality exactly at the separatrix $|\mathbf{L}|^2 = 2TI_2$, where both equal $2TI_2$), and $L_y$ oscillates between $-\beta$ and $+\beta$: the turning value follows directly from conservation by setting $L_x = 0$, which gives $L_y^2 = \beta^2$.

This ODE is the standard form for **Jacobi elliptic-function** parametrisation. The closed-form solution is (Landau-Lifshitz Vol I §37 Eqs. (37.8)–(37.10); Goldstein 2002 §5.6):

$$
\boxed{
\begin{aligned}
L_x(t) \;&=\; \sqrt{\frac{I_1\,(2TI_3 - |\mathbf{L}|^2)}{I_3 - I_1}}\; \mathrm{cn}\big(\Omega^* t + u_0,\, k\big), \\
L_y(t) \;&=\; \sqrt{\frac{I_2\,(2TI_3 - |\mathbf{L}|^2)}{I_3 - I_2}}\; \mathrm{sn}\big(\Omega^* t + u_0,\, k\big), \\
L_z(t) \;&=\; \sqrt{\frac{I_3\,(|\mathbf{L}|^2 - 2TI_1)}{I_3 - I_1}}\; \mathrm{dn}\big(\Omega^* t + u_0,\, k\big),
\end{aligned}
}
$$

with phase $u_0$ fixed by the initial condition. Consistency checks (exact, by the Jacobi identities $\mathrm{cn}^2 = 1 - \mathrm{sn}^2$, $\mathrm{dn}^2 = 1 - k^2\mathrm{sn}^2$): $L_x^2 + L_y^2 + L_z^2 = |\mathbf{L}|^2$ and $L_x^2/I_1 + L_y^2/I_2 + L_z^2/I_3 = 2T$ hold identically for all $t$. The $L_y$ amplitude is $\beta$, matching the turning-point argument above; $L_z$ never vanishes on this branch ($\mathrm{dn} > 0$), which is the polhode-encircles-$\hat{\mathbf{e}}_3$ statement. The characteristic frequency $\Omega^*$ and modulus $k$ are:

$$
\Omega^*
\;=\;
\sqrt{\frac{(I_3 - I_2)\,(|\mathbf{L}|^2 - 2 T I_1)}{I_1 I_2 I_3}},
\qquad
k^2
\;=\;
\frac{(I_2 - I_1)\,(2 T I_3 - |\mathbf{L}|^2)}{(I_3 - I_2)\,(|\mathbf{L}|^2 - 2 T I_1)}.
$$

(The detailed algebra is in Landau-Lifshitz §37 and Goldstein §5.6 with slightly differing conventions; the form here follows Landau-Lifshitz. Implementation note: `scipy.special.ellipk`/`ellipj` take the parameter $m = k^2$, not the modulus $k$.)

> **Correction note (2026-07-11).** The v0.1.0–v0.3.0 published form of this subsection carried erroneous boxed expressions for the reduction prefactor, the amplitude coefficients, $\Omega^*$, and $k^2$, and misidentified the branch (the amplitude of $L_y$ was wrong by a factor ~16 and the implied flip period by 7–63% depending on configuration). The error was found when the flip-period cross-check against a tight-tolerance numerical integration was implemented for the first time (forensic audit 2026-07-11; internal audit record `_audit/forensic-2026-07-11/` in the source monorepo): the corrected forms above match the numerical golden to $\lesssim 10^{-12}$ in period and $8\times10^{-14}$ element-wise in $(L_x, L_y, L_z)(t)$ across the §6.5.5 configurations. Regression guard: `tests/test_elliptic_oracle.py`. This is the fourth boxed-formula intra-doc drift incident (KB-043 class) and the first to have reached a public mirror — precisely because the "cross-check" this box needed was described in fixture metadata but never implemented.

**Period of the elliptic motion.** The trajectory in body-frame $\mathbf{L}$-space is periodic with period:

$$
\boxed{
T_{\mathrm{LL}} \;=\; \frac{4\, K(k)}{\Omega^*},
}
$$

where $K(k) = \int_0^{\pi/2} d\phi / \sqrt{1 - k^2 \sin^2\phi}$ is the complete elliptic integral of the first kind. This is the **closed-form expression for the Dzhanibekov flip period** at finite amplitude — the analytical content of OQ-5.5.

**Dzhanibekov regime: $|\mathbf{L}|^2$ close to $2 T I_2$.** At the unstable intermediate-axis equilibrium ($|\mathbf{L}|^2 \to 2 T I_2$, $k \to 1$), the period diverges: $K(k) \to \infty$ as $k \to 1$ ($K(k) \sim \ln(4/\sqrt{1-k^2})$ as $k \to 1$). Physically: trajectories near the intermediate-axis saddle pass *arbitrarily close* to it on a timescale $\sim (1/\sigma) \ln(1/|\delta\Omega|_{\mathrm{initial}})$, in agreement with the §6.5.2 linearised growth rate (logarithmic dependence on initial perturbation as expected from saddle-point dynamics). For finite but small initial perturbations $|\delta\Omega_x| \ll |\Omega_2^0|$:

$$
T_{\mathrm{flip}} \;\approx\; \frac{2}{\sigma}\, \ln\!\left(\frac{|\Omega_2^0|}{|\delta\Omega|_{\mathrm{initial}}}\right) \;+\; \mathcal{O}(1),
$$

matching the §6.5.2 linearised exponential-growth dynamics in the early-time regime. At large initial perturbations ($|\delta\Omega| \sim |\Omega_2^0|$), the elliptic-modulus $k$ drops well below 1 and the flip period asymptotes to a constant: $T_{\mathrm{flip}} \sim 4 K(k)/\Omega^*$ with $k$ saturating at an inertia-ratio-determined value bounded away from 1. **The flip-period $\to$ amplitude dependence is therefore $\ln(1/\mathrm{amplitude})$ at small amplitudes asymptoting to a constant at large amplitudes** — a falsifiable prediction of the closed-form solution.

**Geometric interpretation.** The Jacobi $\mathrm{sn}, \mathrm{cn}, \mathrm{dn}$ functions parametrise the closed curve at the polhode intersection. As $k \to 0$, they reduce to $\sin, \cos$ and the trajectory becomes the regular precession of an axisymmetric top (cross-check: the $I_1 \to I_2$ limit gives the symmetric top of §6.3, with $\mathrm{sn}, \mathrm{cn}$ → $\sin, \cos$ and modulus $k \to 0$). As $k \to 1$, the period diverges and the trajectory spends most of its time near the saddle — the Dzhanibekov regime. This single closed-form solution interpolates smoothly between the §6.3 free symmetric top and the §6.5.2 Dzhanibekov instability.

### §6.5.4 Podolský anchor (CS-language)

**Partial anchor.** Podolský 2018 covers symmetric tops thoroughly (§8.1 free symmetric top, §8.2 damped, §8.3 Lagrange top) but does not derive the asymmetric-top Jacobi-elliptic solution explicitly. The classical Czech-language reference for this material is Trkal *Mechanika hmotných bodů a tuhého tělesa* (older but standard); modern CS readers may also use Goldstein 2002 in translation if available. The Dzhanibekov effect itself (named after the Soviet cosmonaut Vladimir Dzhanibekov who observed it during his 1985 Salyut-7 mission) is documented in Western space-flight literature (Levi 2014 *Classical Mechanics with Calculus of Variations and Optimal Control* §6 — modern textbook treatment).

### §6.5.5 Numerical tolerance + test-file

| Quantity | Tolerance | Test-file |
|---|---|---|
| Linearised growth rate $\sigma$ vs §6.5.2 analytical formula | $< 5\%$ over the linear regime | `backends/fortran/tests/test_instability_dzhanibekov.f90` (Phase 1 D9 Case 2) + `tests/test_dynamics.py::TestDzhanibekovInstability` (queued engineer-tier) |
| Finite-amplitude flip period vs §6.5.3 boxed $T_{\mathrm{LL}} = 4K(k)/\Omega^*$ | $< 10^{-3}$ relative; degrades systematically near $k \to 1$ saddle approach as $K(k)$ becomes elliptic-modulus-sensitive | `tests/test_elliptic_oracle.py` (added 2026-07-11; the D9 scan validated flips against numerics only, not against this closed form — see §6.5.3 correction note) |
| Logarithmic-amplitude scaling of flip period at small perturbation | $T_{\mathrm{flip}} \cdot \sigma = 2\ln(|\Omega_2^0|/|\delta\Omega|) + \mathcal{O}(1)$ slope match | (same) |
| Conservation $T$, $|\mathbf{L}|^2$ across flips | $< 10^{-11}$ relative as in §6.3 testbed | (same) |
| Number of completed flips over 600 s window | Match the elliptic-function-predicted count, ±1 | (same) |

With §6.5.3's closed-form solution, the full testbed is **Class A** at both linear and finite amplitude. The Jacobi-elliptic solution provides analytical oracles for *every* output of the numerical integration: position-on-polhode, flip period, transition between regimes via $k$-modulus tracking. The previous "Class B numerical-only" designation for the finite-amplitude regime is **retired** — the §6.5.3 closed-form solution resolves the analytical-oracle gap.

## §6.6 Infinite-separation decoupling testbed (Class B)

### §6.6.1 Mathematical setup + expected behaviour

For $|\boldsymbol{\rho}| \to \infty$, the monopole-monopole term $V^{(0)} \to 0$ and both tidal terms $V_{\mathrm{tidal},A}, V_{\mathrm{tidal},B} \to 0$ as $1/|\boldsymbol{\rho}|^3$. The total mutual potential vanishes faster than any inverse-power physical scale. Each body becomes an isolated torque-free top: $\boldsymbol{\tau}_A, \boldsymbol{\tau}_B \to \mathbf{0}$ (proportional to $1/|\boldsymbol{\rho}|^3$), $V \to 0$, and the Lagrangian factors as $L = L_A + L_B$ where each $L_X = T_X$ is the kinetic-only Lagrangian of an isolated rigid body. The system decouples completely.

In practice: place the bodies at large separation $|\boldsymbol{\rho}| = R_{\mathrm{far}}$, choose initial body-frame angular velocities such that each body's rotational dynamics would be independent torque-free Euler precession (Ch 06 §6.3) of the appropriate symmetry class, and verify numerically that the cross-coupling (the residual gravity-gradient torque + back-reaction force at this large but finite $R_{\mathrm{far}}$) stays below the threshold prescribed by the asymptotic scaling.

### §6.6.2 Expected scaling + tolerance

At separation $R_{\mathrm{far}}$, the residual gravity-gradient torque magnitude is $|\boldsymbol{\tau}_A| \sim 3Gm_B \cdot |\mathbf{I}_A - I_0 \mathbb{I}|/R_{\mathrm{far}}^3$. Choosing $R_{\mathrm{far}} = 10^3 \cdot \ell_A$ where $\ell_A$ is the body size, the dimensionless residual coupling is $\sim 10^{-9}$ relative to the torque-free spin angular momentum. Over one orbital period $T = R_{\mathrm{far}}^{3/2}/(GM)^{1/2}$, the accumulated effect on each body's individual angular momentum should stay below $10^{-12}$.

### §6.6.3 Test-file

| Quantity | Tolerance | Test-file |
|---|---|---|
| Cross-coupling between body-$A$ and body-$B$ angular momenta | $< 10^{-12}$ relative per orbital period | `tests/test_dynamics.py::TestInfiniteSeparationDecoupling` (planned, per `_specs/...` §3.2) |
| Each body's individual $|\mathbf{L}_X|$ stability | Class A §6.3 tolerance | (same) |
| Each body's kinetic energy $T_X$ stability | Class A §6.3 tolerance | (same) |

This testbed validates the **multipole-truncation honesty** of Ch 03 §3.4.3 boxed: higher-order multipole terms should decay correctly with $1/r^n$, and no spurious coupling at large $r$ should appear from numerical artifacts.

## §6.7 Cross-implementation bit-identity (Class B procedural — doctrine §4.4 C2)

### §6.7.1 Discipline

When two independent implementations of the boxed Ch 04 §4.7 equation set exist — currently: (a) `dynamics.py` Python first-cut at engineer-tier, (b) `backends/fortran/src/jwst_l2_dynamics.f90` Fortran lege-artis backend at Phase 1 — the cross-implementation outputs must agree at a doctrine-prescribed tolerance. Following the C2 criterion of `_config/LEGE-ARTIS-LANGUAGE-DOCTRINE-v0.1.md` §4.4 (applied at the analytical-mechanics domain), the tolerance is **ULP × √N** per Higham §4.2 Option G, where ULP is the unit in the last place at the magnitude of the compared quantity and $N$ is the number of cumulative arithmetic operations (RK4 steps in our case).

For the first-cut configuration: $N \approx 12{,}000$ RK4 steps over 600 s at $dt = 0.05$ s; ULP at $|\boldsymbol{\rho}| \sim 10$ m is $\sim 10^{-15}\, |\boldsymbol{\rho}|$; ULP × √N $\sim 10^{-13}$ relative. This is the **gate** for the cross-implementation comparison.

### §6.7.2 Testbed implementations

The cross-implementation comparison runs across these Phase 1 deliverables (per `_handoffs/fortran-phase-1/SONNET-BRIEF.md`):

- **D4 — gravity-gradient torque cross-test.** Single configuration; compare Fortran-computed $\boldsymbol{\tau}_X^{\mathrm{body}}$ vs Python-computed at identical inputs. Hard gate at ULP × √N relative tolerance.
- **D8 — symmetric-top trajectory cross-test.** Single 12,000-step RK4 trajectory of the §6.3 testbed configuration; compare full $(t, \boldsymbol{\Omega}_X(t))$ time series. Hard gate at ULP × √N relative tolerance.

After Phase 1 closes (`jwst-l2-fortran-v0.1.0-phase-1-green` tag), the cross-implementation discipline extends to the full equation set with Phase 2 covering the libration + Dzhanibekov tests.

### §6.7.3 What it tests

This is the analytical-mechanics-domain analog of the fourier project's **library-backed vs hand-rolled FFT** discipline (Ch 03 §3.6 + project-owner JWST-L2 captures §3.3): when 2+ implementations agree bit-by-bit, the probability that both share the same algebraic-implementation error decreases dramatically. Cross-implementation bit-identity does not catch errors in the canonical-tier equations themselves (both implementations could implement the same wrong equation); for that, the Class A testbeds against analytical solutions of Ch 06 §6.2–§6.5 are the relevant guard. C2 catches **implementation drift** between the canonical equations and their numerical realisation.

## §6.8 What is established by the end of this chapter

A reader who has worked through §6.1–§6.7 has:

- The **discipline statement** (§6.1): Class A closed-form analytical solutions as oracles for the integrator; Class B edge-case limit configurations as model-domain validators. Both required for honest validation.
- **Class A core testbeds:**
  - §6.2 Kepler limit (analytical Kepler equation in the $\mathbf{I}_X \to \mathbf{0}$ limit; tolerance $10^{-6}$ relative on orbital period).
  - §6.3 Free symmetric top (Podolský §8.1 Eqs. 8.1–8.2; precession-frequency recovery within 0.2 FFT bins; implemented in first-cut with $|dE/E_0| = 6.2 \times 10^{-12}$).
  - §6.4 Gravity-gradient libration linearised (Hughes §8.4 + Podolský §8.3 effective-potential method; libration period within $10^{-4}$ over 10 orbits; Phase 2 testbed).
  - §6.5 Dzhanibekov instability linearised (analytical growth rate $\sigma$; Phase 1 D9 Case 2).
- **Class B core testbeds:**
  - §6.5 Dzhanibekov full-nonlinear (flip dynamics, period amplitude-dependent).
  - §6.6 Infinite-separation decoupling (cross-coupling $< 10^{-12}$ per orbital period).
  - §6.7 Cross-implementation bit-identity (doctrine §4.4 C2 at ULP × √N).
- **Pointer-map** from each symmetry-class specialisation in Ch 05 to its corresponding testbed in this chapter.

**What is NOT yet established at v0.1:**

- The **L2-frame restricted three-body testbeds** (linearised stability eigenvalues at L2; halo orbit closed-form approximation; Lyapunov orbits). v0.2 canonical-tier scope.
- The **symplectic-Lie-group integrator** as an alternative implementation for the C2 cross-implementation gate to validate. The current C2 testbed compares Python ↔ Fortran (two RK4 implementations of the same equations); a symplectic-integrator port would enable a richer C2 cross-validation. v0.2 engineer-tier scope (Ch 04 §4.9 OQ-4.6).
- The **full Lyapunov spectrum** of the asymmetric body in gravity-gradient field (the relevant generalisation of the §6.5 Dzhanibekov testbed to the gravity-gradient-perturbed case, with quantitative Lyapunov-exponent prediction). v0.2 canonical-tier scope (Ch 04 §4.9 OQ-4.8).
- The **relativistic 1PN-corrections** testbed for the regime $v/c \to \mathrm{nonsmall}$. Explicitly out of scope for v0.1 (Ch 00 §6); could become a v0.3+ Class A entry if the scope ever extends.

---

## §6.9 Open questions in this chapter

### Forward-looking (v0.2 scope)

- **OQ-6.1.** Implement and validate the **L2 restricted three-body** testbeds: linearised stability eigenvalues at L2 (Szebehely 1967 §4; Murray & Dermott 1999 §3); halo orbit closed-form approximation (Richardson 1980 third-order expansion; Howell 1984); Lyapunov-orbit family numerical continuation. Each is a Class A or B testbed for the v0.2 canonical tier (which adds the rotating Earth-Sun frame + centrifugal + Coriolis pseudo-forces).
- **OQ-6.2.** Implement the **long-horizon RK4 drift assessment** (`_specs/...` §3.7): integrate §6.3 testbed over $T = 10^4 T_{\mathrm{orbit}}$ and quantify the secular energy drift, confirming the $O(dt^4 \cdot T)$ scaling for the non-symplectic RK4. This is the empirical baseline against which the symplectic-Lie-group integrator (Ch 04 §4.9 OQ-4.6) is compared.
- **OQ-6.3.** Implement the **RK4 order-of-accuracy convergence audit** (`_specs/...` §3.6): integrate §6.3 testbed at $dt \in \{0.5, 0.1, 0.05, 0.01, 0.005, 0.001\}$ s and confirm log-log fit of $|dE|$ vs $dt$ has slope $4.0 \pm 0.1$ in the convergent regime. Doctrine §4.4 C4 analog for ODE integration.
- **OQ-6.4** ***Resolved*** — at §6.4.3 (this chapter): closed-form Lagrange-top + gravity-gradient-libration solution via effective-potential $V_{\mathrm{ef}}(\vartheta)$ method, with cyclic-coordinate first integrals $L_\psi, L_\varphi$, boxed effective potential, closed-form quadrature for $\vartheta(t)$, and three precession regimes (a)/(b)/(c). Falsifiable amplitude-period scaling $T_{\mathrm{lib}}(\vartheta_{\mathrm{amp}}) = (2/\pi)K(k)\cdot T_{\mathrm{lib}}^{\mathrm{linear}}$. Pairs with Ch 05 §5.8 OQ-5.6 resolution.
- **OQ-6.5** ***Resolved*** — at §6.5.3 (this chapter): closed-form Jacobi-elliptic-function parametrisation $(L_x, L_y, L_z) = (\alpha\sqrt{\cdot}\,\mathrm{cn}, \alpha\,\mathrm{sn}, \alpha\sqrt{\cdot}\,\mathrm{dn})$ with elliptic modulus $k$ and characteristic frequency $\Omega^*$ determined by $T$ and $|\mathbf{L}|^2$. Boxed period $T_{\mathrm{LL}} = 4K(k)/\Omega^*$ → Dzhanibekov flip period as falsifiable amplitude-dependent quantity. Pairs with Ch 05 §5.8 OQ-5.5 resolution.

### Procedural / engineer-tier

- **OQ-6.6.** Implement the C2 cross-implementation bit-identity gate (§6.7) as a CI workflow on the GitHub Actions runner for both Python and Fortran backends. Pattern after `lege-artis/fourier`'s v0.2.x cross-backend-bit-identity CI (per fourier doctrine §4.4 C2 testbed-expansion brief).

---

> **Reference anchors.** For the analytical-mechanics testbeds: Landau-Lifshitz Vol I §13–§15 (Kepler), §35 (symmetric top), §37 (asymmetric top); Goldstein 2002 §5.4–§5.7 (rigid-body mechanics including Dzhanibekov + Lagrange top); Hughes 2004 §3–§5 + §8 (spacecraft attitude dynamics + gravity-gradient stabilisation). For the modern Lie-group / variational perspective: Marsden & Ratiu 1999 §13 (Euler-Poincaré); Holm-Schmah-Stoica 2009 §6 ($SO(3)$ rigid body as worked example). For the symplectic / Lie-group integrator framework (Ch 04 §4.9 OQ-4.6 pairing): Hairer-Lubich-Wanner 2006 §VII.5 (symplectic Runge-Kutta on Lie groups). For the L2 restricted three-body extension (v0.2): Szebehely 1967; Murray & Dermott 1999 §3.
>
> **Supplementary anchor for Czech-reading audiences.** Podolský 2018 Kapitola 8 *Aplikace: setrvačníky* provides direct CS-language analytical-limit oracles for three of the six core testbeds: §8.1 (Eqs. 8.1–8.2) for §6.3 free symmetric top; §8.2 damped variant; §8.3 (Eqs. 8.10–8.15) for §6.4 gravity-gradient libration via Lagrange top. The §6.2 Kepler testbed has classical anchors in CS literature outside Podolský (Landau-Lifshitz translation; older Trkal). The §6.5 Dzhanibekov asymmetric-top is *not* directly covered in Podolský — modern CS readers should pair with Goldstein 2002 §5.6.

---

**Next:** v0.2 of the canonical tier (L2 restricted three-body extension; symplectic-Lie-group integrator chapter content; full nonlinear Lagrange-top + Jacobi-elliptic asymmetric-top analytical solutions).

**End 06-analytical-limits-and-validation-gates.md (v0.1, round 3).**
