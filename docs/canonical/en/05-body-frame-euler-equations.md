# Chapter 5 — Body-frame Euler equations as a specialisation

> **Prerequisite reading.** Chapter 02 (kinetic energy on $SE(3) \times SE(3)$, especially §2.7 left-invariance of the rotational metric and §2.7.1 explicit construction). Chapter 03 (gravity-gradient torque, especially §3.7 boxed body-frame expression and §3.7.2 orientation Jacobian + libration-equilibrium stability table after OQ-FORT-1 sign-fix). Chapter 04 (Euler-Lagrange + Euler-Poincaré reduction; this chapter is the *specialisation* of Chapter 04 §4.4 boxed result to the principal-axis frame, with the explicit gravity-gradient torque substituted from Ch 03 §3.7).
>
> **Status.** v0.1 (canonical-tier round 3, 2026-05-25). All four chapter-internal open questions OQ-5.1, OQ-5.2, OQ-5.3, OQ-5.4 resolved in this round; see §5.8 for the close-out summary.

---

## §5.1 What this chapter does

Chapter 04 derived, via Euler-Poincaré reduction of the left-invariant rotational Lagrangian on each $SO(3)$ factor, the general body-frame equations of motion

$$
\mathbf{I}_X \dot{\boldsymbol{\Omega}}_X \;+\; \boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X) \;=\; \boldsymbol{\tau}_X^{\mathrm{body}}, \qquad X \in \{A, B\},
$$

where $\boldsymbol{\Omega}_X$ is the body-frame angular velocity of body $X$, $\mathbf{I}_X$ is its (constant in body frame) inertia tensor, and $\boldsymbol{\tau}_X^{\mathrm{body}}$ is the body-frame external torque on body $X$. **This chapter specialises that general result to the principal-axis frame** (where $\mathbf{I}_X$ is diagonal) and substitutes the explicit gravity-gradient torque from Chapter 03 §3.7 to produce the closed-form scalar equations of motion that the engineer-tier implementation propagates and the Fortran lege-artis backend validates.

The chapter also surfaces the **specific structural features of the principal-axis form** that matter for numerical implementation: the antisymmetric coupling structure of the Euler terms $(I_2 - I_3)\Omega_y\Omega_z$, the special role of the spherical and axisymmetric symmetry classes (where some Euler terms vanish identically), and the connection to the analytical-limit testbeds (free symmetric top, Lagrange top, Dzhanibekov asymmetric-top instability) that Chapter 06 will validate against.

This chapter follows the lege-artis methodology established in Chapters 02–04: every equation is derived from the previous chapter's results explicitly, every symmetry-class specialisation is justified, every analytical limit is named and pointed to its Chapter 06 validation testbed. The chapter is structured to be a **closed specialisation reference** — a reader who knows the result of Chapter 04 (the body-frame Euler equation in coordinate-free form) can come to this chapter for the principal-axis scalar form alone.

## §5.2 Principal-axis specialisation of the body-frame Euler equations

The Ch 04 §4.4 boxed body-frame Euler equation reads, in coordinate-free form:

$$
\mathbf{I}_X \dot{\boldsymbol{\Omega}}_X \;+\; \boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X) \;=\; \boldsymbol{\tau}_X^{\mathrm{body}}.
$$

The inertia tensor $\mathbf{I}_X$ is symmetric and positive-definite (Ch 02 §2.4), so it is diagonalisable by an orthogonal basis change. The **principal-axis frame** is the body-fixed orthonormal frame $(\hat{\mathbf{e}}_1^X, \hat{\mathbf{e}}_2^X, \hat{\mathbf{e}}_3^X)$ in which $\mathbf{I}_X$ takes the diagonal form $\mathbf{I}_X = \mathrm{diag}(I_1^X, I_2^X, I_3^X)$, with the **principal moments of inertia** $I_1^X, I_2^X, I_3^X$ as the diagonal entries (the eigenvalues of $\mathbf{I}_X$). Rotation to the principal-axis frame is a similarity transformation applied to the body-frame equation; we work in that frame throughout this chapter.

### §5.2.1 Component-by-component expansion (resolves OQ-5.1)

Write the boxed equation component-wise in the principal-axis frame. Let $\boldsymbol{\Omega}_X = (\Omega_x, \Omega_y, \Omega_z)$ in the principal-axis basis (suppressing the $X$ subscript on individual components for legibility; $X$ is restored when ambiguous). The two terms on the left are:

**Term 1 — body-frame angular acceleration:**

$$
\mathbf{I}_X \dot{\boldsymbol{\Omega}}_X
\;=\;
\mathrm{diag}(I_1, I_2, I_3)
\begin{pmatrix} \dot\Omega_x \\ \dot\Omega_y \\ \dot\Omega_z \end{pmatrix}
\;=\;
\begin{pmatrix} I_1 \dot\Omega_x \\ I_2 \dot\Omega_y \\ I_3 \dot\Omega_z \end{pmatrix}.
$$

**Term 2 — angular-velocity / angular-momentum cross product:** $\mathbf{I}_X \boldsymbol{\Omega}_X = (I_1 \Omega_x, I_2 \Omega_y, I_3 \Omega_z)$, so:

$$
\boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X)
\;=\;
\begin{vmatrix} \hat{\mathbf{e}}_1 & \hat{\mathbf{e}}_2 & \hat{\mathbf{e}}_3 \\ \Omega_x & \Omega_y & \Omega_z \\ I_1\Omega_x & I_2\Omega_y & I_3\Omega_z \end{vmatrix}
\;=\;
\begin{pmatrix}
\Omega_y(I_3\Omega_z) - \Omega_z(I_2\Omega_y) \\
\Omega_z(I_1\Omega_x) - \Omega_x(I_3\Omega_z) \\
\Omega_x(I_2\Omega_y) - \Omega_y(I_1\Omega_x)
\end{pmatrix}
\;=\;
\begin{pmatrix}
(I_3 - I_2)\,\Omega_y\Omega_z \\
(I_1 - I_3)\,\Omega_z\Omega_x \\
(I_2 - I_1)\,\Omega_x\Omega_y
\end{pmatrix}.
$$

Summing the two terms and equating to the body-frame torque $\boldsymbol{\tau}_X^{\mathrm{body}} = (\tau_x^{\mathrm{body}}, \tau_y^{\mathrm{body}}, \tau_z^{\mathrm{body}})$:

$$
\boxed{
\begin{aligned}
I_1 \dot\Omega_x \;-\; (I_2 - I_3)\, \Omega_y \Omega_z \;&=\; \tau_x^{\mathrm{body}}, \\
I_2 \dot\Omega_y \;-\; (I_3 - I_1)\, \Omega_z \Omega_x \;&=\; \tau_y^{\mathrm{body}}, \\
I_3 \dot\Omega_z \;-\; (I_1 - I_2)\, \Omega_x \Omega_y \;&=\; \tau_z^{\mathrm{body}}.
\end{aligned}
}
$$

These are the classical **Euler dynamic equations of 1758** for a rigid body in its principal-axis frame, recovered here as a structural specialisation of Chapter 04's Euler-Poincaré-derived equation set rather than as a Newton-Euler force-balance statement. The two routes agree (Chapter 04 §4.5 cross-validation); this chapter takes the Euler-Poincaré-derived equation set as the canonical statement and proves that it specialises to the same algebraic form Euler wrote down in 1758. **(Resolves OQ-5.1.)**

The body-$B$ equations are symmetric with $A \leftrightarrow B$ and an independent set of body-frame principal moments and angular velocity components. The two body-frame equation sets are **coupled only through the orientation-dependent gravity-gradient torques** $\boldsymbol{\tau}_X^{\mathrm{body}}$, which depend on the relative orientation of the two bodies through the body-frame direction $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, X}$ (Chapter 03 §3.7 boxed form).

### §5.2.2 Sign-convention sanity check

The cross-product term enters with a **minus** sign on the left-hand side of the boxed equations (equivalently, with a **plus** sign on the right-hand side if it is moved over). This sign comes from the Ch 04 §4.4 derivation of $\mathrm{ad}^*_{\boldsymbol{\Omega}}\,\boldsymbol{\mu} = -\boldsymbol{\Omega} \times \boldsymbol{\mu}$, which propagates through the Euler-Poincaré equation $d(\mathbf{I}\boldsymbol{\Omega})/dt = \mathrm{ad}^*_{\boldsymbol{\Omega}}(\mathbf{I}\boldsymbol{\Omega}) + \boldsymbol{\tau}^{\mathrm{body}}$ to give $\mathbf{I}\dot{\boldsymbol{\Omega}} = -\boldsymbol{\Omega} \times (\mathbf{I}\boldsymbol{\Omega}) + \boldsymbol{\tau}^{\mathrm{body}}$ or equivalently $\mathbf{I}\dot{\boldsymbol{\Omega}} + \boldsymbol{\Omega} \times (\mathbf{I}\boldsymbol{\Omega}) = \boldsymbol{\tau}^{\mathrm{body}}$.

**Cross-check via the Dzhanibekov instability sign.** For an asymmetric body with $I_1 < I_2 < I_3$ rotating principally about the intermediate axis $\hat{\mathbf{e}}_2$ — i.e., $\boldsymbol{\Omega} = (0, \Omega_2^0, 0) + \boldsymbol{\delta\Omega}$ with $|\boldsymbol{\delta\Omega}| \ll |\Omega_2^0|$ — the linearised equations for $\delta\Omega_x, \delta\Omega_z$ derived by substituting into the boxed §5.2.1 equations and keeping first-order terms give an instability eigenvalue $\sigma^2 = (\Omega_2^0)^2 (I_2 - I_1)(I_3 - I_2)/(I_1 I_3) > 0$ (positive, hence unstable). For rotation about either extreme axis ($\hat{\mathbf{e}}_1$ with $\boldsymbol{\Omega} = (\Omega_1^0, 0, 0)$ + perturbation, or $\hat{\mathbf{e}}_3$ similarly), the analogous eigenvalue is $\sigma^2 = -(\Omega_1^0)^2 (I_2 - I_1)(I_3 - I_1)/(I_2 I_3) < 0$ (negative, hence stable). The sign of the boxed Euler equations is therefore consistent with the textbook tennis-racket instability classification (Landau-Lifshitz Vol I §37; Goldstein 2002 §5.6); the Dzhanibekov testbed (Ch 06 §6.5 / `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md` §2.2) is a permanent regression guard.

## §5.3 Symmetry-class specialisations

The structure of the Euler equations simplifies dramatically when the body has additional symmetry. The three relevant cases:

| Body symmetry class | Principal moments | Euler-equation structure |
|---|---|---|
| **Spherical** | $I_1 = I_2 = I_3 \equiv I_0$ | All Euler cross-terms vanish identically; the equations reduce to $I_0 \dot{\boldsymbol{\Omega}} = \boldsymbol{\tau}^{\mathrm{body}}$ — pure torque-driven angular acceleration with no inertial coupling. The torque-free case ($\boldsymbol{\tau} = 0$) gives $\boldsymbol{\Omega} = \mathrm{const}$. Spherical bodies have *no internal rotational dynamics* in the body frame — all complexity is in the orientation evolution. |
| **Axisymmetric** ($I_1 = I_2 \equiv I_\perp$, $I_3 \equiv I_\parallel$) | Two distinct moments | The $\dot\Omega_z$ equation decouples: $I_\parallel \dot\Omega_z = \tau_z^{\mathrm{body}}$. The $(x, y)$ pair couples as $I_\perp \dot\Omega_x = (I_\perp - I_\parallel) \Omega_y \Omega_z + \tau_x^{\mathrm{body}}$ (and cyclic). The torque-free case is the **free symmetric top** with regular precession of $\boldsymbol{\Omega}$ around the symmetry axis at rate $\lambda = (I_\parallel - I_\perp)\Omega_z / I_\perp$. |
| **Asymmetric** (all three $I_i$ distinct) | Three distinct moments | All three Euler cross-terms are nonzero. The torque-free case exhibits the **Dzhanibekov tennis-racket effect**: rotation about the intermediate principal axis ($I_2$ between $I_1$ and $I_3$) is unstable; small perturbations grow into periodic 180° flips of the spin axis. |

### §5.3.1 Spherical case — worked details

In the spherical case $I_1 = I_2 = I_3 = I_0$, every coefficient $(I_i - I_j)$ vanishes, the cross-product term $\boldsymbol{\Omega} \times (\mathbf{I}\boldsymbol{\Omega}) = \boldsymbol{\Omega} \times (I_0 \boldsymbol{\Omega}) = I_0\, \boldsymbol{\Omega} \times \boldsymbol{\Omega} = \mathbf{0}$, and the equation reduces to:

$$
I_0\, \dot{\boldsymbol{\Omega}} \;=\; \boldsymbol{\tau}^{\mathrm{body}}.
$$

This is decoupled across components and identical to Newton's second law applied to a single rotational degree of freedom. In the torque-free case ($\boldsymbol{\tau} = \mathbf{0}$), $\boldsymbol{\Omega}$ is conserved in the body frame, and via $\dot{\mathbf{R}} = \mathbf{R}\,\widehat{\boldsymbol{\Omega}}$ the inertial-frame orientation rotates uniformly about the constant angular-velocity vector. No precession; no nutation; the body is rotationally "boring". In the gravity-gradient case, the gravity-gradient torque vanishes (Ch 03 §3.8.2: for $\mathbf{I}_A = I_0 \mathbb{I}_3$, $\hat{\boldsymbol{\rho}}_{\mathrm{body},A} \times (\mathbf{I}_A \hat{\boldsymbol{\rho}}_{\mathrm{body},A}) = I_0\, \hat{\boldsymbol{\rho}}_{\mathrm{body},A} \times \hat{\boldsymbol{\rho}}_{\mathrm{body},A} = \mathbf{0}$), and the dynamics decouple into pure Kepler orbit + uniform inertial-frame rotation. This is testbed §2.4 / Ch 06 §6.3-limit (with the analytical reduction further reducing to §6.2 Kepler if $B$ is also spherical or a point mass).

### §5.3.2 Axisymmetric case — free symmetric top

For an axisymmetric body $I_1 = I_2 \equiv I_\perp$, $I_3 \equiv I_\parallel$ in the torque-free case, the boxed §5.2.1 equations become:

$$
I_\perp \dot\Omega_x = (I_\perp - I_\parallel)\,\Omega_y \Omega_z, \qquad
I_\perp \dot\Omega_y = -(I_\perp - I_\parallel)\,\Omega_x \Omega_z, \qquad
I_\parallel \dot\Omega_z = 0.
$$

(The $\dot\Omega_z$ equation simplified because $(I_1 - I_2) = 0$ kills the cross-product term.) The third equation says $\Omega_z = \mathrm{const}$ — the body-frame spin about the symmetry axis is conserved. Substituting back into the first two:

$$
\dot\Omega_x = \lambda\, \Omega_y, \qquad \dot\Omega_y = -\lambda\, \Omega_x, \qquad \lambda \equiv \frac{I_\parallel - I_\perp}{I_\perp}\, \Omega_z = \mathrm{const}.
$$

This is the **Euler precession** of the transverse angular-velocity components: $\Omega_x, \Omega_y$ oscillate harmonically at angular frequency $\lambda$. The closed-form solution is:

$$
\Omega_x(t) = A \cos(\lambda t + \phi_0), \qquad \Omega_y(t) = A \sin(\lambda t + \phi_0), \qquad \Omega_z(t) = \mathrm{const},
$$

with $A$ and $\phi_0$ set by initial conditions. The body-frame angular-velocity vector $\boldsymbol{\Omega}(t)$ traces a **body cone** of half-angle $\arctan(A/\Omega_z)$ around the symmetry axis $\hat{\mathbf{e}}_3$ at constant rate $\lambda$. This is the analytical solution that the Phase 1 D8 testbed `backends/fortran/tests/test_symmetric_top.f90` validates numerically, with Podolský 2018 §8.1 Eqs. (8.1)–(8.2) as the direct CS-language analytical oracle.

### §5.3.3 Asymmetric case — Dzhanibekov instability

For an asymmetric body with $I_1 < I_2 < I_3$ in the torque-free case, the boxed §5.2.1 equations have all three cross-product terms nonzero. Three rotational equilibria exist (rotation about each principal axis), with stability classified by:

| Equilibrium | Linearised eigenvalue $\sigma^2$ | Classification |
|---|---|---|
| About $\hat{\mathbf{e}}_1$ ($I_1$, smallest moment) | $-(\Omega_1^0)^2 (I_2 - I_1)(I_3 - I_1)/(I_2 I_3) < 0$ | Stable (oscillatory) |
| About $\hat{\mathbf{e}}_2$ ($I_2$, intermediate moment) | $+(\Omega_2^0)^2 (I_2 - I_1)(I_3 - I_2)/(I_1 I_3) > 0$ | **Unstable** (exponentially growing) |
| About $\hat{\mathbf{e}}_3$ ($I_3$, largest moment) | $-(\Omega_3^0)^2 (I_3 - I_1)(I_3 - I_2)/(I_1 I_2) < 0$ | Stable (oscillatory) |

The unstable middle-axis case is the **Dzhanibekov tennis-racket effect**: small perturbations of an intermediate-axis-aligned spin grow exponentially until amplitude saturates, the spin flips by approximately 180°, then re-stabilises near the opposite intermediate-axis-aligned configuration; periodically thereafter, the body undergoes characteristic flips. The motion is bounded (energy + angular-momentum conservation forces it onto a closed trajectory on the body-frame angular-velocity sphere) but the spin axis flips intermittently. The flip period depends on the initial perturbation amplitude: larger perturbation → faster flips. The complete analytical solution involves Jacobi elliptic functions (Goldstein 2002 §5.6); the linearised instability eigenvalue above gives the e-folding time of the early-time growth.

This is the Phase 1 D9 Case 2 testbed `backends/fortran/tests/test_instability_dzhanibekov.f90`. Each Phase 1 cross-test catches sign discrepancies in the boxed §5.2.1 form by physical-truth-anchored validation (the flip *must* happen, the period *must* match the elliptic-function prediction); a typo in the cross-product term would surface as a wrong-direction flip or as no flip at all.

## §5.4 Substitution of the gravity-gradient torque into the principal-axis form

Insert the body-frame gravity-gradient torque from Chapter 03 §3.7 boxed result,

$$
\boldsymbol{\tau}_X^{\mathrm{body}} \;=\; \frac{3 G m_{\bar X}}{|\boldsymbol{\rho}|^3}\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, X} \times \big(\mathbf{I}_X \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, X}\big),
$$

where $\bar X$ is the *other* body, into the boxed §5.2.1 principal-axis Euler equations. In the principal-axis frame, let $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, X} = (n_1, n_2, n_3)$ (the unit-magnitude inter-body separation direction expressed in body $X$'s principal-axis frame, so $n_1^2 + n_2^2 + n_3^2 = 1$). Then $\mathbf{I}_X \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},X} = (I_1 n_1, I_2 n_2, I_3 n_3)$, and the cross product is:

$$
\hat{\boldsymbol{\rho}}_{\mathrm{body},X} \times (\mathbf{I}_X \hat{\boldsymbol{\rho}}_{\mathrm{body},X})
\;=\;
\begin{vmatrix} \hat{\mathbf{e}}_1 & \hat{\mathbf{e}}_2 & \hat{\mathbf{e}}_3 \\ n_1 & n_2 & n_3 \\ I_1 n_1 & I_2 n_2 & I_3 n_3 \end{vmatrix}
\;=\;
\begin{pmatrix}
(I_3 - I_2)\, n_2 n_3 \\
(I_1 - I_3)\, n_3 n_1 \\
(I_2 - I_1)\, n_1 n_2
\end{pmatrix}.
$$

The gravity-gradient torque in principal-axis components is therefore:

$$
\boxed{
\boldsymbol{\tau}_X^{\mathrm{body}}
\;=\;
\frac{3 G m_{\bar X}}{|\boldsymbol{\rho}|^3}\,
\begin{pmatrix}
(I_3^X - I_2^X)\, n_2 n_3 \\
(I_1^X - I_3^X)\, n_3 n_1 \\
(I_2^X - I_1^X)\, n_1 n_2
\end{pmatrix}.
}
$$

The right-hand side carries the **same $(I_i - I_j)$-type prefactors** as the Euler cross-terms on the left-hand side of the boxed §5.2.1 equations. This is not a coincidence: both come from the antisymmetric / commutator structure of $\mathfrak{so}(3)$ acting on a symmetric $\mathbf{I}$, and both vanish identically for a spherical body (all $I_i - I_j = 0$).

**Combined closed-form principal-axis equations of motion** (substituting the boxed gravity-gradient torque into the boxed §5.2.1 form):

$$
\begin{aligned}
I_1 \dot\Omega_x - (I_2 - I_3)\,\Omega_y\Omega_z \;&=\; (3Gm_{\bar X}/|\boldsymbol{\rho}|^3)\,(I_3 - I_2)\, n_2 n_3, \\
I_2 \dot\Omega_y - (I_3 - I_1)\,\Omega_z\Omega_x \;&=\; (3Gm_{\bar X}/|\boldsymbol{\rho}|^3)\,(I_1 - I_3)\, n_3 n_1, \\
I_3 \dot\Omega_z - (I_1 - I_2)\,\Omega_x\Omega_y \;&=\; (3Gm_{\bar X}/|\boldsymbol{\rho}|^3)\,(I_2 - I_1)\, n_1 n_2.
\end{aligned}
$$

This is the **closed-form principal-axis equation set** the engineer-tier `dynamics.py` Euler-equations branch (Ch 03 §3.9 table → `dynamics.py::state_derivative` rotational lines) propagates, and the Fortran lege-artis backend `backends/fortran/src/jwst_l2_dynamics.f90` implements at Phase 1.

### §5.4.1 Libration-equilibrium specialisation (resolves OQ-5.2 via OQ-FORT-1 cross-check)

Specialise to an axisymmetric body $A$ ($I_1 = I_2 = I_\perp$, $I_3 = I_\parallel$, with principal axes aligned with body frame) at the **radial-pointing equilibrium**: the body's symmetry axis aligned with $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$, so $\hat{\boldsymbol{\rho}}_{\mathrm{body},A} = (0, 0, 1)$ → $(n_1, n_2, n_3) = (0, 0, 1)$ at equilibrium. Substituting into the boxed §5.4 gravity-gradient torque:

$$
\boldsymbol{\tau}_A^{\mathrm{body}}\big|_{\mathrm{eq}}
\;=\;
\frac{3Gm_B}{|\boldsymbol{\rho}|^3}\,
\big((I_3 - I_2)\cdot 0 \cdot 1, \; (I_1 - I_3)\cdot 1 \cdot 0, \; (I_2 - I_1)\cdot 0 \cdot 0\big)
\;=\;
\mathbf{0}.
$$

The equilibrium is a true critical point: zero net torque. **Now perturb the orientation by a small angle $\theta$ about, say, body-axis $\hat{\mathbf{e}}_1$.** Adopt the Ch 03 §3.7 convention for the body-frame infinitesimal rotation: $\mathbf{R}_q^A \to \mathbf{R}_q^A(\mathbb{I} + \widehat{\delta\boldsymbol{\theta}_A})$ with $\delta\boldsymbol{\theta}_A = (\theta, 0, 0)$. Then by the Ch 03 §3.7 identity $\delta\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = -\delta\boldsymbol{\theta}_A \times \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$:

$$
\delta\hat{\boldsymbol{\rho}}_{\mathrm{body},A}
\;=\;
-(\theta, 0, 0) \times (0, 0, 1)
\;=\;
-\theta\,(\hat{\mathbf{e}}_1 \times \hat{\mathbf{e}}_3)
\;=\;
-\theta\,(-\hat{\mathbf{e}}_2)
\;=\;
+\theta\,\hat{\mathbf{e}}_2.
$$

So the perturbed body-frame separation direction is $\hat{\boldsymbol{\rho}}_{\mathrm{body},A}(\theta) \approx (0, +\theta, 1)$, giving $(n_1, n_2, n_3) = (0, +\theta, 1)$, and the perturbed torque is:

$$
\boldsymbol{\tau}_A^{\mathrm{body}}(\theta)
\;\approx\;
\frac{3Gm_B}{|\boldsymbol{\rho}|^3}\,
\big((I_3 - I_2)\cdot (+\theta) \cdot 1, \; 0, \; 0\big)
\;=\;
\frac{3Gm_B(I_\parallel - I_\perp)}{|\boldsymbol{\rho}|^3}\,(\theta, 0, 0),
$$

using $I_2 = I_\perp, I_3 = I_\parallel$ for the axisymmetric body. The libration equation for the small-angle perturbation is then $I_\perp\, \ddot\theta = \tau_x^{\mathrm{body}}(\theta)$:

$$
\boxed{
I_\perp\, \ddot\theta \;=\; \frac{3Gm_B(I_\parallel - I_\perp)}{|\boldsymbol{\rho}|^3}\, \theta.
}
$$

This is a linear ODE $\ddot\theta = k\, \theta$ with $k = 3Gm_B(I_\parallel - I_\perp)/(I_\perp\,|\boldsymbol{\rho}|^3)$. The sign of $k$ determines stability:

- **Prolate** ($I_\parallel < I_\perp$, $k < 0$): $\ddot\theta = -|k|\theta$ → simple harmonic oscillator → **stable libration** at frequency $\omega_{\mathrm{lib}}^2 = |k| = 3Gm_B(I_\perp - I_\parallel)/(I_\perp\,|\boldsymbol{\rho}|^3) > 0$. The symmetry-axis-along-radial attitude is the gravity-gradient-stabilised attitude (the body "wants" to point its long axis at the primary).
- **Oblate** ($I_\parallel > I_\perp$, $k > 0$): $\ddot\theta = +k\theta$ → exponentially growing → **unstable** with e-folding time $\tau_{\mathrm{inst}} = 1/\sqrt{k}$. The stable attitude has the symmetry axis *perpendicular* to $\hat{\boldsymbol{\rho}}$ (a flat sunshield-like body aligns face-on).

**Comparison with Ch 03 §3.7.2 orientation Jacobian.** Ch 03 §3.7.2 boxed (after OQ-FORT-1 sign-fix 2026-05-24 PM) gives the libration-equilibrium spring-constant matrix:

$$
\mathbf{J}^{\mathrm{orient}}_{A}\big|_{\mathrm{eq}}
\;=\; \alpha\, (I_\parallel - I_\perp)\, \mathrm{diag}(1, 1, 0)
\;=\; \frac{3Gm_B(I_\parallel - I_\perp)}{|\boldsymbol{\rho}|^3}\,\mathrm{diag}(1,1,0),
$$

with the eigenvalue $J_{11} = 3Gm_B(I_\parallel - I_\perp)/|\boldsymbol{\rho}|^3$. Multiplying the linearised relation $\delta\boldsymbol{\tau}_A^{\mathrm{body}} = \mathbf{J}^{\mathrm{orient}}_A\, \delta\boldsymbol{\theta}_A$ into the libration equation $I_\perp\ddot\theta = \delta\tau_x = J_{11}\theta$:

$$
I_\perp\,\ddot\theta \;=\; J_{11}\, \theta \;=\; \frac{3Gm_B(I_\parallel - I_\perp)}{|\boldsymbol{\rho}|^3}\,\theta,
$$

which is **identical** to this section's boxed result. The two routes — (a) substituting the boxed §5.4 principal-axis-form torque into the boxed Ch 05 §5.2.1 Euler equation at the libration equilibrium (this §5.4.1), and (b) applying the Ch 03 §3.7.2 orientation Jacobian at the equilibrium directly — produce the same linearised libration equation with the same stability classification (prolate stable / oblate unstable). **(Resolves OQ-5.2.)**

This is the **boxed-formula sanity-check substitution discipline** (per KB-OPUS-DOC-INTRA-INCONSISTENCY): the §5.4 boxed principal-axis-form gravity-gradient torque is verified against the independently-derived Ch 03 §3.7.2 boxed orientation Jacobian by substituting a specific worked case (axisymmetric body, on-axis equilibrium, transverse perturbation) and confirming agreement. The Phase 1 D5 testbed `backends/fortran/tests/test_orientation_jacobian.f90` (T2a libration-equilibrium eigenvalues) is the cross-implementation regression guard at ULP × √N tolerance — same data, two derivation routes, three implementations (canonical-tier Ch 03, canonical-tier Ch 05, engineer-tier Fortran D5) all in agreement.

### §5.4.2 Symmetry-class structure of the substituted equations

Reading the symmetry-class table of §5.3 in conjunction with the boxed §5.4 gravity-gradient torque:

| Body symmetry class | Number of nonzero principal-axis torque components | Dynamical content |
|---|---|---|
| Spherical ($I_1 = I_2 = I_3$) | 0 (all $(I_i - I_j) = 0$) | Gravity-gradient torque vanishes identically; rotation decouples from orbit; testbed §2.4 / Ch 06 §6.3-limit. |
| Axisymmetric ($I_1 = I_2 \neq I_3$) | 2 generically (the $\tau_z^{\mathrm{body}}$ vanishes when $I_1 = I_2$, leaving $\tau_x, \tau_y$ as the libration-drive components) | Two-dimensional libration plane perpendicular to the symmetry axis; one cyclic direction along the axis. Testbed §2.5 / Ch 06 §6.4 (small-angle libration / Lagrange top). |
| Asymmetric (all $I_i$ distinct) | 3 generically | Three independent libration / instability modes; the gravity-gradient torque couples to all three orientation degrees of freedom. Multi-mode coupling; non-closed-form except in special cases. Testbed §2.5 + §3.5 / Ch 06 §6.4-extended. |

The structure is **mirror-symmetric** with the symmetry-class table of §5.3 for the torque-free case: the same anisotropy condition $(I_i - I_j) \neq 0$ that makes the inertia cross-term active in the torque-free Euler equation makes the gravity-gradient torque component active in the substituted equation. This is no accident: both terms originate in the antisymmetric tensor structure of $\boldsymbol{\Omega} \times (\mathbf{I}\boldsymbol{\Omega})$ vs $\hat{\boldsymbol{\rho}}_{\mathrm{body}} \times (\mathbf{I}\hat{\boldsymbol{\rho}}_{\mathrm{body}})$, and they have parallel rank tables.

## §5.5 Constants of motion and structural conservation laws

The body-frame Euler equations have well-known constants of motion in the torque-free case:

- **Kinetic energy in the body frame** $T = \tfrac{1}{2}(I_1 \Omega_x^2 + I_2 \Omega_y^2 + I_3 \Omega_z^2)$ is conserved.
- **Magnitude of body-frame angular momentum** $|\mathbf{L}|^2 = (I_1 \Omega_x)^2 + (I_2 \Omega_y)^2 + (I_3 \Omega_z)^2$ is conserved.

**Verification.** Differentiate $T$ with respect to time:

$$
\frac{dT}{dt} \;=\; I_1 \Omega_x \dot\Omega_x + I_2 \Omega_y \dot\Omega_y + I_3 \Omega_z \dot\Omega_z.
$$

Substituting the torque-free boxed §5.2.1 Euler equations ($I_i \dot\Omega_i = (I_j - I_k)\Omega_j\Omega_k$ for the appropriate cyclic indices):

$$
\frac{dT}{dt}
\;=\;
\Omega_x \cdot (I_2 - I_3)\Omega_y\Omega_z + \Omega_y \cdot (I_3 - I_1)\Omega_z\Omega_x + \Omega_z \cdot (I_1 - I_2)\Omega_x\Omega_y
\;=\;
\Omega_x \Omega_y \Omega_z \cdot \big[(I_2 - I_3) + (I_3 - I_1) + (I_1 - I_2)\big]
\;=\;
0.
$$

The three coefficients sum to zero identically, confirming $T = \mathrm{const}$. The analogous calculation for $|\mathbf{L}|^2 = (I_1\Omega_x)^2 + (I_2\Omega_y)^2 + (I_3\Omega_z)^2$ proceeds:

$$
\frac{d|\mathbf{L}|^2}{dt} = 2\, I_1\Omega_x \cdot I_1 \dot\Omega_x + \ldots = 2\, I_1\Omega_x \cdot (I_2 - I_3)\Omega_y\Omega_z + \ldots = 2\Omega_x\Omega_y\Omega_z\, [I_1(I_2 - I_3) + I_2(I_3 - I_1) + I_3(I_1 - I_2)] = 0,
$$

again by the algebraic identity $I_1(I_2 - I_3) + I_2(I_3 - I_1) + I_3(I_1 - I_2) = 0$ (expand and cancel).

### §5.5.1 Geometric interpretation

The torque-free trajectory in body-frame $\boldsymbol{\Omega}$-space is the **intersection** of two ellipsoids:

- The **kinetic-energy ellipsoid** $\tfrac{1}{2}(I_1 \Omega_x^2 + I_2 \Omega_y^2 + I_3 \Omega_z^2) = T_0$ (a fixed energy surface).
- The **angular-momentum-magnitude ellipsoid** $(I_1\Omega_x)^2 + (I_2\Omega_y)^2 + (I_3\Omega_z)^2 = L_0^2$ (a fixed angular-momentum-squared surface).

Both are smooth ellipsoids in $\boldsymbol{\Omega}$-space. Their intersection is a 1-dimensional curve (generically; degenerate cases at the special symmetry-class fixed points). The trajectory $\boldsymbol{\Omega}(t)$ traces this curve; it is closed (because both surfaces are bounded) and periodic except at the unstable intermediate-axis fixed points of §5.3.3 where the curve passes through a saddle. The closed-form solution involves Jacobi elliptic functions; the analytical detail is in Landau-Lifshitz Vol I §37 and Goldstein 2002 §5.6. **(Resolves OQ-5.3 partially — the torque-free constants are established; the externally-torqued case is now treated.)**

### §5.5.2 With external torque — total-energy conservation

When the body experiences external torque (our gravity-gradient case), kinetic energy $T$ is **not** conserved alone. The work-energy theorem gives:

$$
\frac{dT}{dt} \;=\; \boldsymbol{\Omega} \cdot \boldsymbol{\tau}^{\mathrm{body}},
$$

which is the power delivered to the rotational degree of freedom by the external torque. (Derivation: differentiate $T = \tfrac{1}{2}\boldsymbol{\Omega}^T \mathbf{I} \boldsymbol{\Omega}$; substitute the body-frame Euler equation $\mathbf{I}\dot{\boldsymbol{\Omega}} = -\boldsymbol{\Omega}\times(\mathbf{I}\boldsymbol{\Omega}) + \boldsymbol{\tau}^{\mathrm{body}}$; use $\boldsymbol{\Omega}\cdot[\boldsymbol{\Omega}\times(\mathbf{I}\boldsymbol{\Omega})] = 0$ scalar-triple-product cyclic identity.)

The torque is the body-frame projection of $-\partial V/\partial \boldsymbol{\theta}$; the corresponding **potential-energy time derivative** is:

$$
\frac{dV_{\mathrm{tidal},X}}{dt}
\;=\;
\frac{\partial V_{\mathrm{tidal},X}}{\partial q_X} \cdot \dot q_X
\;+\; \frac{\partial V_{\mathrm{tidal},X}}{\partial \boldsymbol{\rho}} \cdot \dot{\boldsymbol{\rho}}
\;=\;
-\boldsymbol{\tau}_X^{\mathrm{body}} \cdot \boldsymbol{\Omega}_X
\;+\; \nabla_{\boldsymbol{\rho}} V_{\mathrm{tidal},X}\cdot \dot{\boldsymbol{\rho}},
$$

with the first term being the rate at which the torque does negative work on the body (the body gains $T$ at the rate the orientation potential loses $V$), and the second the rate at which the back-reaction force (Ch 04 §4.6.2 boxed) does work on the relative orbit. Summing:

$$
\frac{d(T + V)}{dt} \;=\; 0 \quad \text{(total energy conservation)}.
$$

The full conservation law is therefore **total mechanical energy** $E = T + V$, with $V$ given by Ch 03 §3.4.3 boxed truncated mutual potential. The Ch 04 §4.6.4 Newton's-third-law balance plus this rotational work-energy relation give:

- $\boldsymbol{P}_{\mathrm{CoM}}$ conserved (translational invariance — Ch 02 §2.6.1 + Noether).
- $\boldsymbol{L}_{\mathrm{total}}$ conserved (rotational invariance — Ch 02 + Ch 04 §4.6.4).
- $E = T + V$ conserved (time-translation invariance + the work-energy / potential-energy bookkeeping above).

These are the **three structural conservation laws** any honest integration of the boxed §4.7 equation set must respect. The first-cut engineer-tier diagnostics report $|dE/E_0| = 6.2 \times 10^{-12}$ and $|d\mathbf{L}/L_0| = 2.1 \times 10^{-10}$ over 600 seconds at $dt = 0.05$ s — the small but nonzero drift is the cumulative RK4 integration error (non-symplectic; $O(dt^4 \cdot T)$ scaling), not a violation of the structural conservation laws. **(Resolves OQ-5.3 fully — torque-free constants $(T, |\mathbf{L}|^2)$ and torqued total-energy conservation $E = T + V$ both established with explicit work-energy bookkeeping.)**

## §5.6 Cross-references to validation testbeds (Chapter 06 pointers)

Each symmetry-class specialisation of §5.3 has an analytical-limit solution that becomes a validation testbed in Chapter 06:

- **Spherical body, torque-free** → trivial validation: $\boldsymbol{\Omega} = \mathrm{const}$, orientation rotates uniformly. Ch 06 §6.3-spherical-limit / `_specs/...` §3.4 zero-rotation initial conditions extended case.
- **Axisymmetric body, torque-free** → **free symmetric top** with regular precession; testbed validates the body-frame solution against the analytical formula $\Omega_x(t) = A \cos(\lambda t + \phi_0)$, $\Omega_y(t) = A \sin(\lambda t + \phi_0)$, $\Omega_z = \mathrm{const}$. This is the Phase 1 D8 test `backends/fortran/tests/test_symmetric_top.f90`, anchored to Podolský §8.1 Eqs. (8.1)–(8.2) as the CS-language analytical oracle. Ch 06 §6.3.
- **Asymmetric body, torque-free** → **Dzhanibekov tennis-racket effect**; testbed validates the predicted finite-amplitude periodic flips. Phase 1 D9 Case 2 test `backends/fortran/tests/test_instability_dzhanibekov.f90`. Ch 06 §6.5.
- **Axisymmetric body, gravity-gradient torque (small libration)** → linearised libration spring constant matches Chapter 03 §3.7.2 result $\alpha(I_\parallel - I_\perp)$; full-amplitude trajectory matches the **Lagrange-top** effective-potential analysis. Phase 2 testbed; anchored to Podolský §8.3 Eqs. (8.10)–(8.15) as the CS-language analytical oracle. Ch 06 §6.4.
- **Asymmetric body, gravity-gradient torque** → multi-mode libration with stability classified by the full Jacobian eigenstructure (Chapter 03 §3.7.2). Phase 2 testbed; no closed-form analytical solution, validates via cross-implementation bit-identity (doctrine §4.4 C2). Ch 06 §6.5-extended.

Chapter 06 brings each testbed to a quantitative pass/fail gate.

## §5.7 What is established by the end of this chapter

A reader who has worked through §5.1–§5.6 has:

- The classical Euler dynamic equations of 1758 (boxed §5.2.1) derived by component-by-component expansion of Chapter 04 §4.4's coordinate-free Euler-Poincaré-derived body-frame Lagrangian — confirming the structural specialisation. **OQ-5.1 resolved.**
- The sign-convention sanity check (§5.2.2) via the Dzhanibekov instability eigenvalue, consistent with textbook tennis-racket-theorem classification.
- Explicit **symmetry-class specialisations** (§5.3): spherical (no internal rotational dynamics), axisymmetric (free symmetric top with closed-form precession solution at body-cone half-angle $\arctan(A/\Omega_z)$ and precession rate $\lambda = (I_\parallel - I_\perp)\Omega_z/I_\perp$), and asymmetric (Dzhanibekov instability about the intermediate principal axis with explicit eigenvalue formulas at the three principal-axis equilibria).
- The **gravity-gradient torque substituted explicitly** into the principal-axis form (§5.4 boxed); structural mirror symmetry observed between the inertia-anisotropy structure of the torque-free Euler cross-terms and that of the gravity-gradient torque components.
- The **libration-equilibrium specialisation** (§5.4.1) reproducing the Ch 03 §3.7.2 boxed orientation Jacobian eigenvalue $J_{11} = 3Gm_B(I_\parallel - I_\perp)/|\boldsymbol{\rho}|^3$ by an independent derivation route. **OQ-5.2 resolved** via the boxed-formula-sanity-check-substitution discipline; libration analysis from this chapter's substituted equations agrees with Ch 03 §3.7.2 orientation Jacobian analysis. Cross-implementation regression guard at Phase 1 D5 testbed.
- The **constants-of-motion structure** (§5.5.1 + §5.5.2): in the torque-free case, body-frame kinetic energy and angular-momentum-squared are conserved by the identity $I_1(I_2-I_3) + I_2(I_3-I_1) + I_3(I_1-I_2) = 0$; geometric interpretation as intersection of energy-ellipsoid and angular-momentum-magnitude-ellipsoid in body-frame $\boldsymbol{\Omega}$-space. In the externally-torqued case, total mechanical energy $E = T + V$ is conserved, with explicit work-energy / potential-energy bookkeeping. **OQ-5.3 resolved.**
- The **pointer-map** from each symmetry class to the corresponding Chapter 06 analytical-limit validation testbed (§5.6).

**What is NOT yet established:**

- The complete validation against analytical-limit testbeds at quantitative pass/fail gates. **Chapter 06.**
- The symplectic-Lie-group integrator that would replace RK4 for long-horizon accurate integration of these equations. **v0.2 engineer-tier scope** (Ch 04 §4.9 OQ-4.6).
- The L2 / restricted-three-body extension adding centrifugal and Coriolis pseudo-forces to the gravity-gradient torque expression. **v0.2 canonical tier.**
- The full nonlinear stability analysis at the libration equilibria via the energy-momentum method. **v0.2 canonical-tier scope** (Ch 04 §4.9 OQ-4.7).

---

## §5.8 Open questions in this chapter

### Resolved this round (canonical-tier round 3, 2026-05-25)

- **OQ-5.1** ***Resolved*** — at §5.2.1: explicit component-by-component expansion of $\boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X)$ in the principal-axis frame produces the boxed scalar Euler equations $I_1 \dot\Omega_x - (I_2 - I_3)\Omega_y\Omega_z = \tau_x^{\mathrm{body}}$ (and cyclic). Recovered from Ch 04 §4.4 coordinate-free result by structural specialisation.
- **OQ-5.2** ***Resolved*** — at §5.4.1: libration-equilibrium specialisation produces libration frequency $\omega^2 = 3Gm_B(I_\parallel - I_\perp)/(I_\perp |\boldsymbol{\rho}|^3)$ via direct linearisation of the boxed §5.4 substituted Euler equations. Spring-constant magnitude matches Ch 03 §3.7.2 boxed orientation Jacobian eigenvalue $J_{11} = \alpha(I_\parallel - I_\perp)$ exactly. Sign convention (prolate stable / oblate destabilising) consistent with Ch 03 §3.7.2 stability table after OQ-FORT-1 sign-fix. Boxed-formula sanity-check substitution discipline confirmed via the agreement between Ch 03 and Ch 05 independent derivations.
- **OQ-5.3** ***Resolved*** — at §5.5.1 + §5.5.2: torque-free constants $T$ and $|\mathbf{L}|^2$ conserved by explicit algebraic identity (sum of $(I_i - I_j)$ cycles equals zero); torqued-case total energy $E = T + V$ conserved via work-energy / potential-energy bookkeeping, with explicit derivation at §5.5.2. Phase 1 conservation diagnostics $|dE/E_0|, |d\mathbf{L}/L_0|$ are the engineer-tier regression guard.
- **OQ-5.4** ***Resolved*** — at §5.6 + cross-check against Ch 04 §4.7.1 table: the boxed §5.4 substituted principal-axis equations are line-by-line consistent with `dynamics.py::state_derivative` rotational lines (162–163) + the Phase 1 Fortran `backends/fortran/src/jwst_l2_dynamics.f90` Euler-equation propagation. The cross-tier match is exact at the rotational level; the only gap is the back-reaction force on the orbit (§5.4 substituted equations describe rotational dynamics only; the back-reaction is Ch 04 §4.6.2 + §4.7.1 / §4.9 OQ-4.5 territory, not Ch 05).

### Resolved in this round's extension (canonical-tier round 3+, 2026-05-25)

- **OQ-5.5** ***Resolved*** — at Ch 06 §6.5.3: closed-form Jacobi-elliptic-function solution of the torque-free asymmetric-top problem authored. Boxed parametrisation $(L_x, L_y, L_z) = (\alpha\sqrt{\cdot}\,\mathrm{cn}(\Omega^*t, k), \alpha\,\mathrm{sn}(\Omega^*t, k), \alpha\sqrt{\cdot}\,\mathrm{dn}(\Omega^*t, k))$ with characteristic frequency $\Omega^*$ and elliptic modulus $k$ derived from $T$ and $|\mathbf{L}|^2$. Boxed period $T_{\mathrm{LL}} = 4K(k)/\Omega^*$ → Dzhanibekov flip period as a falsifiable function of initial perturbation amplitude. Smooth interpolation between §6.3 free symmetric top ($k \to 0$) and the unstable intermediate-axis equilibrium ($k \to 1$). Anchor: Landau-Lifshitz Vol I §37, Goldstein 2002 §5.6.
- **OQ-5.6** ***Resolved*** — at Ch 06 §6.4.3: closed-form Lagrange-top + gravity-gradient-libration solution via the effective-potential $V_{\mathrm{ef}}(\vartheta)$ method. Boxed effective potential $V_{\mathrm{ef}}(\vartheta) = (L_\varphi - L_\psi\cos\vartheta)^2/(2I_\perp\sin^2\vartheta) - 3Gm_B(I_\perp - I_\parallel)\cos^2\vartheta/(2|\boldsymbol{\rho}|^3)$. Cyclic-coordinate first integrals $L_\psi, L_\varphi$ derived; closed-form quadrature for $\vartheta(t)$ given; three precession regimes (a)/(b)/(c) classified per Podolský §8.3 phenomenology. Finite-amplitude libration period $T_{\mathrm{lib}}(\vartheta_{\mathrm{amp}}) = (2/\pi)K(k)\cdot T_{\mathrm{lib}}^{\mathrm{linear}}$ as elliptic-function-of-amplitude correction. Anchor: Podolský 2018 §8.3 Eqs. (8.10)–(8.15) directly; Goldstein 2002 §5.7.

---

> **Reference anchor.** Euler 1758 — *Théorie du mouvement des corps solides ou rigides* (St. Petersburg Academy) — the original derivation of the body-frame equations bearing his name. Modern presentations: Goldstein 2002 §5.4 (textbook treatment with principal-axis discussion); Landau-Lifshitz 1976 §35 (concise tensor-form derivation); Arnold 1989 §29 (modern Lie-group perspective, connecting to Chapter 04's Euler-Poincaré framing); Hughes 2004 §3.3 (spacecraft-attitude domain treatment, with explicit symmetry-class specialisations and the principal-axis-aligned-with-radial libration analysis that anchors §5.6).
>
> **Supplementary anchor for Czech-reading audiences.** Podolský 2018 Kapitola 7 §7.2 *Eulerovy dynamické rovnice* is the direct CS-language analog of this chapter's central boxed result (§5.2.1): his **Eq. (7.15) is byte-identical to our §5.2.1 boxed Euler equations**, derived in his presentation via the Newton-Euler $d\mathbf{L}/dt = \mathbf{M}$ route plus the transport equation (6.8). The two routes (Newton-Euler in Podolský §7.2; Euler-Poincaré reduction in our Chapter 04) produce the same equations — same agreement-as-validation pattern as in our §3.6 and §4.5. For the validation-testbed connections (§5.6), Podolský Kapitola 8 *Aplikace: setrvačníky* is the direct CS-language reference: §8.1 free symmetric top + §8.2 setrvačník se třením (damped) + §8.3 Lagrange top with $V_{\mathrm{ef}}$ method. Each Phase 1 / Phase 2 test in our backend ports has its CS-language analytical oracle in Podolský Kap 8.

---

**Next chapter:** [06 — Analytical limits and validation gates](06-analytical-limits-and-validation-gates.md)

**End 05-body-frame-euler-equations.md (v0.1, round 3).**
