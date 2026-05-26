# Chapter 2 — Kinetic energy on $SE(3) \times SE(3)$ in intrinsic form

> **Prerequisite reading.** Chapter 01 (configuration manifold and variational principle). Familiarity with the inertia tensor as a multilinear object is helpful but is rebuilt here from scratch for completeness.
>
> **Status.** v0.1.

---

## §2.1 What this chapter does

Chapter 1 set the configuration manifold $Q = SE(3) \times SE(3)$ and identified the tangent bundle $TQ$ on which the Lagrangian lives. This chapter writes down the **kinetic energy** half of that Lagrangian explicitly, derived from first principles by integration over the mass distribution of each rigid body. The result is a smooth function $T: TQ \to \mathbb{R}$ that is quadratic in velocities and, viewed geometrically, defines a **Riemannian metric** on $TQ$. The potential energy half (Chapter 03) will close the construction; assembly of the full Lagrangian and derivation of the equations of motion follow in Chapter 04.

The derivation is **intrinsic**: kinetic energy is a physical scalar, the same in any frame. We express it in coordinates only when computation demands it (the body frame turns out to be the simplest coordinate choice for the rotational part, the inertial frame for the translational part). This intrinsic / coordinate-expression distinction matters because it survives all subsequent reductions and specialisations.

## §2.2 Kinetic energy of one rigid body — from the mass distribution

Consider a single rigid body $A$ with mass distribution $\rho_A(\mathbf{x})$ supported on some compact region $\mathcal{B}_A \subset \mathbb{R}^3$. The body's configuration in the inertial frame is described (per Chapter 1, §1.2) by a centre-of-mass position $\mathbf{R}_A \in \mathbb{R}^3$ and an orientation $q_A \in SO(3)$. A material point of the body labelled by its body-frame position $\mathbf{x} \in \mathcal{B}_A$ (taken with respect to the body's centre of mass) sits in the inertial frame at:

$$
\mathbf{r}_A(\mathbf{x}, t) = \mathbf{R}_A(t) + \mathbf{R}_q^A(t)\, \mathbf{x},
$$

where $\mathbf{R}_q^A$ is the rotation matrix associated with $q_A$. Differentiating with respect to time, the inertial velocity of that material point is:

$$
\mathbf{v}_A(\mathbf{x}, t) = \dot{\mathbf{R}}_A + \dot{\mathbf{R}}_q^A\, \mathbf{x} = \dot{\mathbf{R}}_A + \hat{\boldsymbol{\omega}}_A\, \mathbf{R}_q^A\, \mathbf{x},
$$

using the relation $\dot{\mathbf{R}}_q^A = \hat{\boldsymbol{\omega}}_A\, \mathbf{R}_q^A$ from Chapter 1, §1.4, where $\boldsymbol{\omega}_A \in \mathbb{R}^3$ is the angular velocity of body $A$ expressed in the inertial frame.

The kinetic energy of body $A$ is the volume integral of $\tfrac{1}{2} \rho |\mathbf{v}|^2$ over the body:

$$
T_A = \frac{1}{2} \int_{\mathcal{B}_A} \rho_A(\mathbf{x})\, |\mathbf{v}_A(\mathbf{x}, t)|^2\, dV.
$$

Substituting the velocity expression and expanding the square:

$$
|\mathbf{v}_A|^2 = |\dot{\mathbf{R}}_A|^2 + 2\, \dot{\mathbf{R}}_A \cdot \big(\hat{\boldsymbol{\omega}}_A\, \mathbf{R}_q^A\, \mathbf{x}\big) + |\hat{\boldsymbol{\omega}}_A\, \mathbf{R}_q^A\, \mathbf{x}|^2.
$$

Integrate term by term. The first term gives $|\dot{\mathbf{R}}_A|^2 \int_{\mathcal{B}_A} \rho_A\, dV = m_A\, |\dot{\mathbf{R}}_A|^2$, where $m_A$ is the total mass. The second term contains the integral $\int \rho_A\, \mathbf{x}\, dV$, which **vanishes identically** because $\mathbf{x}$ is measured from the centre of mass: by definition $\int \rho_A\, \mathbf{x}\, dV = \mathbf{0}$. This is the algebraic content of the **centre-of-mass theorem** (König's theorem).

The third term requires care. Using the antisymmetry of $\hat{\boldsymbol{\omega}}$ and the identity $|\mathbf{a} \times \mathbf{b}|^2 = |\mathbf{a}|^2 |\mathbf{b}|^2 - (\mathbf{a} \cdot \mathbf{b})^2$, with $\hat{\boldsymbol{\omega}}\, \mathbf{y} = \boldsymbol{\omega} \times \mathbf{y}$:

$$
|\hat{\boldsymbol{\omega}}_A\, \mathbf{R}_q^A\, \mathbf{x}|^2 = |\boldsymbol{\omega}_A \times (\mathbf{R}_q^A\, \mathbf{x})|^2 = |\boldsymbol{\omega}_A|^2 |\mathbf{R}_q^A\, \mathbf{x}|^2 - (\boldsymbol{\omega}_A \cdot \mathbf{R}_q^A\, \mathbf{x})^2.
$$

Since $\mathbf{R}_q^A$ is orthogonal, $|\mathbf{R}_q^A\, \mathbf{x}| = |\mathbf{x}|$ and $\boldsymbol{\omega}_A \cdot \mathbf{R}_q^A\, \mathbf{x} = (\mathbf{R}_q^A)^T \boldsymbol{\omega}_A \cdot \mathbf{x} = \boldsymbol{\Omega}_A \cdot \mathbf{x}$ where $\boldsymbol{\Omega}_A = (\mathbf{R}_q^A)^T \boldsymbol{\omega}_A$ is the angular velocity in the **body frame** (introduced in Chapter 1, §1.4). Hence:

$$
|\hat{\boldsymbol{\omega}}_A\, \mathbf{R}_q^A\, \mathbf{x}|^2 = |\boldsymbol{\omega}_A|^2 |\mathbf{x}|^2 - (\boldsymbol{\Omega}_A \cdot \mathbf{x})^2.
$$

Carrying out the integration, the third term becomes:

$$
\int_{\mathcal{B}_A} \rho_A \left[ |\boldsymbol{\omega}_A|^2 |\mathbf{x}|^2 - (\boldsymbol{\Omega}_A \cdot \mathbf{x})^2 \right] dV.
$$

This is recognisably $\boldsymbol{\Omega}_A^T\, \mathbf{I}_A\, \boldsymbol{\Omega}_A$ where $\mathbf{I}_A$ is the **inertia tensor of body $A$** in its body frame (defined precisely in §2.4 below); using $|\boldsymbol{\omega}_A| = |\boldsymbol{\Omega}_A|$ (orthogonal transformations preserve norm), the identity holds in either frame.

## §2.3 König's decomposition

Collecting terms from §2.2, the kinetic energy of body $A$ decomposes as:

$$
\boxed{
T_A \;=\; \underbrace{\tfrac{1}{2}\, m_A\, |\dot{\mathbf{R}}_A|^2}_{\textstyle T_A^{\mathrm{trans}}}
\;+\; \underbrace{\tfrac{1}{2}\, \boldsymbol{\Omega}_A^T\, \mathbf{I}_A\, \boldsymbol{\Omega}_A}_{\textstyle T_A^{\mathrm{rot}}}.
}
$$

This is **König's theorem** for the kinetic energy of a rigid body: total kinetic energy equals the kinetic energy of the centre-of-mass motion plus the rotational kinetic energy about the centre of mass. The cross-term vanished because of the centre-of-mass choice for $\mathbf{x}$.

Two consequences worth flagging immediately:

1. The translational and rotational parts **decouple at the level of the kinetic energy**. Any coupling between translation and rotation must come from the potential energy (Chapter 03) — specifically from the gravity-gradient term, which depends on both position and orientation.

2. The rotational part is naturally expressed in the **body frame** ($\boldsymbol{\Omega}_A$, $\mathbf{I}_A$), not the inertial frame. The inertia tensor in the body frame is time-independent (the body is rigid; its mass distribution does not change relative to its own frame). The same quantity in the inertial frame, $\mathbf{R}_q^A\, \mathbf{I}_A (\mathbf{R}_q^A)^T$, is time-dependent through the rotation matrix. This asymmetry is why body-frame coordinates are preferred for rotational dynamics — it is not a convention, it is the only choice that keeps the inertia tensor static.

> **Reference anchor.** Landau-Lifshitz Vol I §32 derives König's theorem in a closely parallel form (their notation: $\mathbf{V}$ for centre-of-mass velocity, $\boldsymbol{\Omega}$ for body-frame angular velocity; same content).

## §2.4 The inertia tensor

The **inertia tensor** $\mathbf{I}_A$ of body $A$ in its body frame is the symmetric positive-definite $3 \times 3$ matrix:

$$
(\mathbf{I}_A)_{ij} \;=\; \int_{\mathcal{B}_A} \rho_A(\mathbf{x}) \left(|\mathbf{x}|^2\, \delta_{ij} - x_i\, x_j\right) dV.
$$

The diagonal entries are the moments of inertia about each body-frame axis; the off-diagonal entries are the products of inertia. By the spectral theorem for symmetric matrices, $\mathbf{I}_A$ admits an orthonormal eigenbasis: the **principal axes** of the body. In the principal-axis body frame, $\mathbf{I}_A = \mathrm{diag}(I_{xx,A}, I_{yy,A}, I_{zz,A})$ is diagonal with non-negative entries.

**Parallel-axis theorem.** If a body is composed of subcomponents $\{C_k\}$ with inertia tensors $\mathbf{I}_{C_k}$ about their own centres of mass at body-frame positions $\mathbf{c}_k$, the composite body's inertia tensor about its centre of mass is:

$$
\mathbf{I}_A \;=\; \sum_k \left[ \mathbf{I}_{C_k} + m_{C_k} \left(|\mathbf{c}_k|^2\, \mathbb{I}_3 - \mathbf{c}_k \mathbf{c}_k^T\right) \right].
$$

The bracketed second term is the parallel-axis contribution: it is the inertia tensor of a point mass $m_{C_k}$ at position $\mathbf{c}_k$ relative to the composite-body centre of mass.

> **Code cross-reference.** `geometries.py::parallel_axis(I_centroid, mass, r_offset)` implements exactly this formula using a $\mathbf{r}\mathbf{r}^T$ outer product. The composite-body inertia primitives `disk_inertia`, `cylinder_inertia`, `thin_rod_inertia`, `solid_cone_inertia` provide $\mathbf{I}_{C_k}$ values for standard subcomponent shapes. The first-cut `geometries.py` exposes three reference bodies used throughout the canonical-tier worked examples: (i) **`make_jwst_like()`** — simplified JWST-style 4-component composite (sunshield + bus + boom + primary), total mass 2450 kg, produces $\mathbf{I} = \mathrm{diag}(2.33 \times 10^4,\ 2.33 \times 10^4,\ 1.54 \times 10^4)\ \text{kg m}^2$ (max/min anisotropy 1.52, **prolate** since component z-spread drives the parallel-axis contribution to $I_\perp$ above the local $I_\parallel$); used as the C2 cross-test oracle for the composite-construction tests. (ii) **`make_probe_like()`** — small 2-component cylindrical-plus-dish body, total mass 100 kg, $\mathbf{I} = \mathrm{diag}(29.06,\ 29.06,\ 5.00)\ \text{kg m}^2$ (anisotropy 5.81, prolate). (iii) **`make_oblate_reference_body()`** — synthetic single-cylinder body sized by construction to produce exactly $\mathbf{I} = \mathrm{diag}(2540,\ 2540,\ 3870)\ \text{kg m}^2$ (anisotropy 1.52, **oblate**); used as the reference body for the gravity-gradient instability worked example in §3.7.2 and the dI/dt sensitivity example in §2.5.1 below. See `_handoffs/fortran-phase-1/INVESTIGATION-OQ-FORT-6.md` for the history of why two anisotropy-1.52 bodies (one prolate, one oblate) coexist in this project.

The eigenvalues of $\mathbf{I}_A$ (the principal moments) are real and non-negative. The **triangle inequality of moments** $I_{xx} + I_{yy} \ge I_{zz}$ (and cyclic) is a structural fact about mass distributions in 3D space; any "inertia tensor" that violates it does not correspond to any physical mass distribution.

> **Reference anchor.** Goldstein 2002 §5.3 (inertia tensor as a Cartesian tensor of rank 2); Landau-Lifshitz Vol I §32 (principal axes and moments); Hughes 2004 Ch 4 (parallel-axis theorem for spacecraft, including composite-body computation patterns essentially identical to the `geometries.py` approach). Thorne & Blandford 2017 §1.1.1 ("The Geometric Viewpoint on the Laws of Physics") and §1.3 ("Tensor Algebra without a Coordinate System") frame the inertia tensor in coordinate-free language as a vector-valued linear function $\mathbf{I}(\,\cdot\,, \boldsymbol{\omega}) = \mathbf{J}$ (angular momentum as the image of angular velocity under $\mathbf{I}$); this is the geometric-viewpoint complement to the coordinate-component derivation given above. T&B does *not* derive the inertia tensor from a mass-distribution integral in this volume — the rigid-body-mechanics material belongs to a companion volume not bundled with the uploaded edition (Modern Classical Physics: Optics, Fluids, Plasmas, Elasticity, Relativity, and Statistical Physics).

## §2.5 Body-frame vs. inertial-frame angular velocity: bookkeeping

The relationship $\boldsymbol{\Omega} = \mathbf{R}_q^T \boldsymbol{\omega}$ between body-frame and inertial-frame angular velocity introduced in Chapter 1, §1.4 has a few subtleties worth pinning down before Chapter 3.

In the body frame, the inertia tensor $\mathbf{I}_A$ is **constant in time** (the body is rigid). The rotational kinetic energy is:

$$
T_A^{\mathrm{rot}} = \tfrac{1}{2}\, \boldsymbol{\Omega}_A^T\, \mathbf{I}_A\, \boldsymbol{\Omega}_A.
$$

In the inertial frame, the inertia tensor at time $t$ is $\mathbf{I}_A^{\mathrm{inertial}}(t) = \mathbf{R}_q^A(t)\, \mathbf{I}_A\, (\mathbf{R}_q^A(t))^T$, which is **time-dependent through $\mathbf{R}_q^A(t)$**. The rotational kinetic energy is then:

$$
T_A^{\mathrm{rot}} = \tfrac{1}{2}\, \boldsymbol{\omega}_A^T\, \mathbf{I}_A^{\mathrm{inertial}}(t)\, \boldsymbol{\omega}_A.
$$

Algebraically, these are the same number (substitute $\boldsymbol{\omega}_A = \mathbf{R}_q^A\, \boldsymbol{\Omega}_A$ and the time-dependence of $\mathbf{I}_A^{\mathrm{inertial}}$ cancels against the rotation of $\boldsymbol{\omega}_A$).

The **angular momentum** has the analogous bookkeeping:

$$
\mathbf{L}_A^{\mathrm{body}} = \mathbf{I}_A\, \boldsymbol{\Omega}_A, \qquad
\mathbf{L}_A^{\mathrm{inertial}} = \mathbf{R}_q^A\, \mathbf{L}_A^{\mathrm{body}} = \mathbf{I}_A^{\mathrm{inertial}}(t)\, \boldsymbol{\omega}_A.
$$

The **inertial-frame** angular momentum is the conserved quantity when external torques are zero. The body-frame angular momentum is generally not conserved; it satisfies Euler's equations (derived in Chapter 5). This conservation asymmetry is a direct consequence of which frame is inertial.

> **Code cross-reference.** `dynamics.py::inertia_inertial_frame(I_body, R_q)` rotates the body-frame inertia tensor to the inertial frame via $\mathbf{R}\, \mathbf{I}\, \mathbf{R}^T$ similarity. `dynamics.py::total_angular_momentum` computes the inertial-frame $\mathbf{L}_A$ used for the conservation diagnostic, then sums across bodies. The first-cut achieves $|d\mathbf{L}/\mathbf{L}_0| = 2.1 \times 10^{-10}$ over 600 s, the round-off floor for double-precision arithmetic at this integration horizon.

### §2.5.1 Time evolution of $\mathbf{I}_A^{\mathrm{inertial}}(t)$: off-diagonal sensitivity (engineer-tier bridge)

The body-frame inertia tensor $\mathbf{I}_A$ is constant; the inertial-frame inertia tensor $\mathbf{I}_A^{\mathrm{inertial}}(t) = \mathbf{R}_q^A(t)\, \mathbf{I}_A\, (\mathbf{R}_q^A(t))^T$ is time-dependent through the orientation. Landau-Lifshitz and the textbook tradition treat this as obvious and move on. For implementation of real calculations — output representation, sensitivity analysis, integrator step-size selection — the explicit **time derivative** $d\mathbf{I}_A^{\mathrm{inertial}}/dt$ and the behaviour of its **off-diagonal elements** are the natural bridge from the canonical formalism to engineering practice. The discussion that follows is therefore the engineer-tier complement to §2.5, requested in project-owner direction 2026-05-23.

Differentiating $\mathbf{I}_A^{\mathrm{inertial}}(t)$ in time and using $\dot{\mathbf{R}}_q^A = \widehat{\boldsymbol{\omega}_A^{\mathrm{inertial}}}\, \mathbf{R}_q^A$ from Chapter 1, §1.4 (with $\boldsymbol{\omega}_A^{\mathrm{inertial}} = \mathbf{R}_q^A \boldsymbol{\Omega}_A$ the inertial-frame angular velocity of body $A$):

$$
\frac{d\mathbf{I}_A^{\mathrm{inertial}}}{dt}
\;=\;
\dot{\mathbf{R}}_q^A\, \mathbf{I}_A\, (\mathbf{R}_q^A)^T
\;+\;
\mathbf{R}_q^A\, \mathbf{I}_A\, (\dot{\mathbf{R}}_q^A)^T
\;=\;
\widehat{\boldsymbol{\omega}_A^{\mathrm{inertial}}}\, \mathbf{I}_A^{\mathrm{inertial}}
\;-\;
\mathbf{I}_A^{\mathrm{inertial}}\, \widehat{\boldsymbol{\omega}_A^{\mathrm{inertial}}},
$$

i.e., the matrix commutator:

$$
\boxed{
\frac{d\mathbf{I}_A^{\mathrm{inertial}}}{dt}
\;=\;
\big[\,\widehat{\boldsymbol{\omega}_A^{\mathrm{inertial}}}\, ,\, \mathbf{I}_A^{\mathrm{inertial}}\,\big].
}
$$

This is the **rotational transport equation** for the symmetric tensor $\mathbf{I}_A^{\mathrm{inertial}}$. It is the analog, for a tensor, of the elementary $d\mathbf{v}/dt = \boldsymbol{\omega} \times \mathbf{v}$ for a body-fixed vector. The trace of $\mathbf{I}_A^{\mathrm{inertial}}$ is preserved (commutators are traceless), as is the eigenvalue spectrum (similarity transformations are eigenvalue-preserving) — consistent with the principal moments being intrinsic body properties.

**Component form, principal-axis frame at the evaluation instant.** Consider an instant at which $\mathbf{R}_q^A = \mathbb{I}$, so $\mathbf{I}_A^{\mathrm{inertial}} = \mathbf{I}_A = \mathrm{diag}(I_1, I_2, I_3)$ is diagonal in inertial-frame coordinates. With body-frame angular velocity $\boldsymbol{\omega} = (\omega_1, \omega_2, \omega_3)$ (equal to the inertial-frame angular velocity at this instant), direct computation of $[\widehat{\boldsymbol{\omega}}, \mathbf{I}_A]_{ij}$ gives the off-diagonal-element rates as:

$$
\left.\frac{d (\mathbf{I}_A^{\mathrm{inertial}})_{ij}}{dt}\right|_{\mathbf{R}_q^A = \mathbb{I}}
\;=\;
\varepsilon_{ijk}\, \omega_k\, (I_i - I_j), \qquad i \neq j.
$$

The diagonal-element rates vanish at this instant: $(d\mathbf{I}/dt)_{ii} = 0$ (which must be true if the principal moments are conserved). The formula is antisymmetric in $(i,j)$, consistent with the symmetry of $\mathbf{I}_A^{\mathrm{inertial}}$ at all times.

**Sensitivity, by body symmetry class.** The off-diagonal rate formula reveals immediately how the body's symmetries propagate into the inertial-frame inertia tensor's time evolution:

| Body symmetry class | Principal moments | Off-diagonal rates at $\mathbf{R}_q^A = \mathbb{I}$ |
|---|---|---|
| Spherical | $I_1 = I_2 = I_3 \equiv I_0$ | All zero. $\mathbf{I}^{\mathrm{inertial}} = I_0 \mathbb{I}_3$ is invariant under rotation; no orientation sensitivity in any quantity that couples to $\mathbf{I}^{\mathrm{inertial}}$. |
| Axisymmetric ($I_1 = I_2 \neq I_3$) | One distinct moment | $(1,2)$ rate vanishes; $(1,3)$ and $(2,3)$ rates nonzero only when $\boldsymbol{\omega}$ has perpendicular components $\omega_2$ or $\omega_1$ respectively. |
| Asymmetric (all distinct) | Three distinct moments | All six off-diagonal rates can be simultaneously nonzero. |

**Numerical sanity check, oblate reference body.** Use the synthetic oblate axisymmetric reference body `make_oblate_reference_body()` (see Chapter 2 §2.4 code cross-reference + OQ-FORT-6 resolution 2026-05-24): $\mathbf{I}_A = \mathrm{diag}(2540, 2540, 3870)\ \text{kg m}^2$ (oblate, with $I_\parallel = I_3 = 3870 > I_\perp = I_1 = I_2 = 2540$, so $I_\parallel / I_\perp = 1.52$). For a slow-tumbling spacecraft with $\boldsymbol{\omega} = (0.01, 0, 0)\ \text{rad/s}$ (rotation period $T_{\mathrm{rot}} = 2\pi/|\boldsymbol{\omega}| = 628\ \text{s}$), the formula above predicts:
- $(d\mathbf{I}^{\mathrm{inertial}}/dt)_{2,3} = \omega_1 (I_3 - I_2) = 0.01 \cdot 1330 = 13.3\ \text{kg m}^2 / \text{s}$,

with the other off-diagonal rates zero by the axisymmetric structure. The time for this off-diagonal element to reach a typical diagonal magnitude $\sim 2540\ \text{kg m}^2$ is $\tau_{\mathrm{off}} \sim 2540 / 13.3 \approx 190\ \text{s} \approx T_{\mathrm{rot}}/3.3$ — within the first-cut integration window of $600\ \text{s}$, the off-diagonal element traces out approximately one full oscillation cycle, comparable in amplitude to the diagonal elements. The off-diagonal elements **are** the body's orientation, represented in inertial-frame components.

> **Two reference bodies, on purpose.** The canonical-tier work uses two complementary axisymmetric reference bodies. (i) `make_oblate_reference_body()` — synthetic single-cylinder body with $\mathbf{I} = \mathrm{diag}(2540, 2540, 3870)$ by construction, **oblate** ($I_\parallel > I_\perp$); used here and in Chapter 03 §3.7.2 for worked-example numerics. (ii) `make_jwst_like()` — simplified 4-component JWST-shaped composite (sunshield + bus + boom + primary mirror) with $\mathbf{I} = \mathrm{diag}(2.33 \times 10^4, 2.33 \times 10^4, 1.54 \times 10^4)$, **prolate** ($I_\parallel < I_\perp$) due to the parallel-axis contribution from the components' z-spread; used for composite-construction demonstrations and as the C2 cross-test oracle. The qualitative comparison between the two is a natural cross-validation: at the gravity-gradient libration equilibrium, the oblate body is unstable (§3.7.2) and the prolate body is stable (§3.8.3). See `_handoffs/fortran-phase-1/INVESTIGATION-OQ-FORT-6.md` for the history.

**Bridge to gravity-gradient torque sensitivity.** The body-frame gravity-gradient torque derived in Chapter 03 §3.7,
$\boldsymbol{\tau}_A^{\mathrm{body}} = (3 G m_B / |\boldsymbol{\rho}|^3)\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A})$, depends only on body-frame quantities (constant $\mathbf{I}_A$ and orientation-dependent $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$). The inertial-frame view, equivalent and useful for output representation, expresses the same torque as
$\boldsymbol{\tau}_A^{\mathrm{inertial}} = (3 G m_B / |\boldsymbol{\rho}|^3)\, \hat{\boldsymbol{\rho}} \times (\mathbf{I}_A^{\mathrm{inertial}}(t) \cdot \hat{\boldsymbol{\rho}})$ with $\hat{\boldsymbol{\rho}}$ in inertial frame. In this form the torque's time evolution and its **sensitivity to initial conditions** are visible directly:
- **Sensitivity to initial mutual position** $\boldsymbol{\rho}(0)$: linearised, $\delta \boldsymbol{\tau}_A \propto \delta \boldsymbol{\rho} \cdot \nabla \boldsymbol{\tau}_A$, and the Jacobian $\nabla \boldsymbol{\tau}_A$ involves $\mathbf{I}_A^{\mathrm{inertial}}(0)$ — the off-diagonals of which encode the body's orientation relative to the inertial axes.
- **Sensitivity to body orientation** $q_A(0)$: small rotations $\delta q_A$ rotate $\mathbf{I}_A^{\mathrm{inertial}}$ by $\delta \mathbf{I}^{\mathrm{inertial}} = [\widehat{\delta q_A}, \mathbf{I}^{\mathrm{inertial}}]$ — the same commutator structure as $d\mathbf{I}/dt$. The torque response is first-order in this off-diagonal perturbation.
- **Two-body assembly sensitivity.** The total dynamics couples through both bodies' orientation-time-evolution and their mutual $\boldsymbol{\rho}$. For two near-spherical bodies, the dynamics is approximately decoupled (gravity-gradient torque vanishes in the spherical limit; §3.8.2 reduction). For two highly anisotropic bodies in close mutual orientation alignment, the assembly's effective anisotropy multiplies and the dynamics is correspondingly more sensitive to initial conditions — the relevant figure of merit is the **product of the two anisotropies** $(I_\parallel/I_\perp)_A \cdot (I_\parallel/I_\perp)_B$ weighted by $\hat{\boldsymbol{\rho}}$-alignment.

**Engineer-tier integrator implication.** The off-diagonal time-scale $\tau_{\mathrm{off}}$ derived above sets a lower bound on the resolved time-scale of the rotational dynamics. For RK4 to track these oscillations accurately, the integration step $dt$ must satisfy $dt \ll \tau_{\mathrm{off}}$. The first-cut configuration uses $dt = 0.05\ \text{s}$ against $\tau_{\mathrm{off}} \approx 190\ \text{s}$, giving roughly four orders of margin — comfortable. For an extreme-anisotropy body or a very rapid rotation, the margin shrinks and the step-size requirement tightens; this is the kind of practical-calculation constraint that the engineer-tier bridge surfaces and that the canonical formalism alone does not make visible.

> **Reference anchor.** The rotational transport equation $d\mathbf{T}/dt = [\widehat{\boldsymbol{\omega}}, \mathbf{T}]$ for a symmetric tensor in a rotating frame is standard; see Hughes 2004 §4 (spacecraft attitude dynamics, formulated for the inertia tensor specifically) and Goldstein 2002 §4.7-§4.10 (rates of change in rotating frames, general tensor case). The off-diagonal sensitivity analysis given here in terms of body-symmetry classes is implicit in Hughes Ch 8 (gravity-gradient stabilisation) but is more naturally read by the engineer-tier audience when written explicitly.

> **Supplementary anchor for Czech-reading audiences.** Podolský 2018 *Teoretická mechanika ve třech knihách* (Matfyzpress, ISBN 978-80-7378-499-7) Kapitola 6 §6.3 derives the same transport equation in coordinate-free form, with the boxed result Eq. (6.8): $\left(\frac{d\mathbf{w}}{dt}\right)\!\big|_{\mathrm{prostor}} = \left(\frac{d\mathbf{w}}{dt}\right)\!\big|_{\mathrm{t\check{e}leso}} + \boldsymbol{\Omega} \times \mathbf{w}$ for any vector $\mathbf{w}$. Applied column-wise to the columns of $\mathbf{I}^{\mathrm{inertial}}$ (each of which is a body-fixed vector expressed in inertial coordinates), the vector form (6.8) yields the matrix-commutator form boxed above. The auxiliary algebra at Podolský Eqs. (6.2) $\boldsymbol{\Omega}' = \frac{dA}{dt} A^t$ (body-frame) and (6.10) $\boldsymbol{\Omega} = A^t \frac{dA}{dt}$ (spatial-frame) provides the explicit trivialisation that makes the body-vs-inertial frame distinction algebraically transparent. Native-language anchor for the CS canonical-tier translation (`docs/canonical/cs/`, forthcoming).

## §2.6 Total kinetic energy of the two-body system

The two-body kinetic energy is the sum:

$$
T \;=\; T_A + T_B
\;=\; \tfrac{1}{2}\, m_A\, |\dot{\mathbf{R}}_A|^2 + \tfrac{1}{2}\, m_B\, |\dot{\mathbf{R}}_B|^2
+ \tfrac{1}{2}\, \boldsymbol{\Omega}_A^T\, \mathbf{I}_A\, \boldsymbol{\Omega}_A
+ \tfrac{1}{2}\, \boldsymbol{\Omega}_B^T\, \mathbf{I}_B\, \boldsymbol{\Omega}_B.
$$

There are **no cross-terms** between body $A$ and body $B$ in the kinetic energy. Each body contributes its own translational + rotational kinetic energy, independent of the other. Geometrically: kinetic energy on the product manifold $TQ = T(SE(3)_A \times SE(3)_B) = TSE(3)_A \times TSE(3)_B$ decomposes as a sum of kinetic energies on each factor.

This is **not a special property** of the two-body system; it is a general fact about kinetic energy of mechanical systems. Kinetic energy is a sum over particles (or, equivalently, over rigid-body subsystems treated as composite particles). The notion of "mechanical coupling" between two systems is always a statement about their potential energy, never their kinetic energy, in classical mechanics.

(In Lagrangian mechanics with **velocity-dependent potentials** — e.g., a charged particle in a magnetic field where $L = \tfrac{1}{2} m \dot{\mathbf{x}}^2 + q \mathbf{A}(\mathbf{x}) \cdot \dot{\mathbf{x}}$ — the linear-in-velocity term mimics a kinetic cross-term, but it is structurally a potential-energy contribution. For our purely gravitational two-body system, no such term arises and the kinetic energy is the pure quadratic-in-velocity sum above.)

We can rewrite the translational part using **centre-of-mass + relative coordinates** (the convention adopted in Chapter 1, §1.3.1). Let $\mathbf{R}_{\mathrm{CoM}} = (m_A \mathbf{R}_A + m_B \mathbf{R}_B)/(m_A + m_B)$ and $\boldsymbol{\rho} = \mathbf{R}_B - \mathbf{R}_A$. The translational kinetic energy becomes:

$$
\tfrac{1}{2}\, m_A\, |\dot{\mathbf{R}}_A|^2 + \tfrac{1}{2}\, m_B\, |\dot{\mathbf{R}}_B|^2
= \tfrac{1}{2}\, M\, |\dot{\mathbf{R}}_{\mathrm{CoM}}|^2 + \tfrac{1}{2}\, \mu\, |\dot{\boldsymbol{\rho}}|^2,
$$

where $M = m_A + m_B$ is the total mass and $\mu = m_A m_B / M$ is the **reduced mass**. In the centre-of-mass frame, $\dot{\mathbf{R}}_{\mathrm{CoM}} = \mathbf{0}$ identically (this is the frame choice locked in Chapter 1) and the translational kinetic energy reduces to the relative-motion term $\tfrac{1}{2}\, \mu\, |\dot{\boldsymbol{\rho}}|^2$.

In the centre-of-mass frame, the total kinetic energy is therefore:

$$
\boxed{
T_{\mathrm{CoM}} \;=\; \tfrac{1}{2}\, \mu\, |\dot{\boldsymbol{\rho}}|^2
\;+\; \tfrac{1}{2}\, \boldsymbol{\Omega}_A^T\, \mathbf{I}_A\, \boldsymbol{\Omega}_A
\;+\; \tfrac{1}{2}\, \boldsymbol{\Omega}_B^T\, \mathbf{I}_B\, \boldsymbol{\Omega}_B.
}
$$

This is the form of the kinetic energy that enters the Lagrangian $L = T - V$ in Chapter 04. The change of variables $(\mathbf{R}_A, \mathbf{R}_B) \to (\mathbf{R}_{\mathrm{CoM}}, \boldsymbol{\rho})$ has been used informally above; the next subsection makes its canonical-mechanics status explicit.

### §2.6.1 The centre-of-mass decomposition as a symplectic point transformation

The change of variables above is more than a notational convenience: it is a **canonical (symplectic) point transformation** on the configuration manifold's cotangent bundle. Making this explicit is what resolves OQ-2.1 of the working-spec OQ list — it pins down the formal status of a coordinate change that the kinetic-energy derivation uses informally.

**Configuration-space transformation.** Define the linear map on $\mathbb{R}^3 \times \mathbb{R}^3$:

$$
\begin{pmatrix} \mathbf{R}_{\mathrm{CoM}} \\ \boldsymbol{\rho} \end{pmatrix}
\;=\;
\mathbf{J}\, \begin{pmatrix} \mathbf{R}_A \\ \mathbf{R}_B \end{pmatrix},
\qquad
\mathbf{J} \;=\;
\begin{pmatrix} m_A/M & m_B/M \\ -\mathbb{I}_3 & +\mathbb{I}_3 \end{pmatrix}
\,\otimes\, \mathbb{I}_3
\;\equiv\;
\begin{pmatrix} m_A/M & m_B/M \\ -1 & +1 \end{pmatrix},
$$

where the second equality drops the trivial $\otimes \mathbb{I}_3$ factor and the matrix entries are block scalars acting independently on each spatial component. By direct computation, $\det \mathbf{J} = m_A/M + m_B/M = 1$. The transformation is therefore volume-preserving on the configuration manifold.

**Cotangent lift to $T^*Q$.** A standard theorem in symplectic geometry (Marsden & Ratiu 1999 §6.3; Arnold 1989 §3.16) states that any diffeomorphism $\Phi: Q \to Q'$ extends uniquely to a symplectic diffeomorphism $T^*\Phi: T^*Q \to T^*Q'$ — the **cotangent lift** — under which the canonical 1-form $\theta = p_i\, dq^i$ is preserved. For a linear point transformation $Q' = \mathbf{J}\, Q$, the lift acts on momenta by the inverse-transpose Jacobian: $\mathbf{p}' = (\mathbf{J}^T)^{-1} \mathbf{p}$.

Computing $\mathbf{J}^{-1}$ and its transpose explicitly:

$$
\mathbf{J}^{-1} \;=\; \begin{pmatrix} 1 & -m_B/M \\ 1 & +m_A/M \end{pmatrix},
\qquad
(\mathbf{J}^T)^{-1} = (\mathbf{J}^{-1})^T \;=\; \begin{pmatrix} 1 & 1 \\ -m_B/M & m_A/M \end{pmatrix}.
$$

The induced momentum transformation $(\mathbf{p}_A, \mathbf{p}_B) \to (\mathbf{P}_{\mathrm{CoM}}, \mathbf{p}_{\mathrm{rel}})$ is therefore:

$$
\boxed{
\mathbf{P}_{\mathrm{CoM}} \;=\; \mathbf{p}_A + \mathbf{p}_B, \qquad
\mathbf{p}_{\mathrm{rel}} \;=\; \frac{-m_B \mathbf{p}_A + m_A \mathbf{p}_B}{M}.
}
$$

**Verification: $\mathbf{p}_{\mathrm{rel}} = \mu \dot{\boldsymbol{\rho}}$ as expected.** Substituting $\mathbf{p}_A = m_A \dot{\mathbf{R}}_A$, $\mathbf{p}_B = m_B \dot{\mathbf{R}}_B$:

$$
\mathbf{p}_{\mathrm{rel}}
\;=\; \frac{-m_B (m_A \dot{\mathbf{R}}_A) + m_A (m_B \dot{\mathbf{R}}_B)}{M}
\;=\; \frac{m_A m_B}{M}\, (\dot{\mathbf{R}}_B - \dot{\mathbf{R}}_A)
\;=\; \mu\, \dot{\boldsymbol{\rho}}.
$$

The conjugate momentum to the relative-position coordinate is reduced-mass times relative velocity, as required.

**Symplecticity check.** The cotangent-lift theorem above guarantees that $\theta$ is preserved, hence so is $\omega = -d\theta$. A direct check makes this concrete: substitute the inverse transformations $\mathbf{p}_A = (\mathbf{p}_{\mathrm{rel}}\, m_A - \mathbf{P}_{\mathrm{CoM}}\, m_A)/(-M) = \mathbf{P}_{\mathrm{CoM}}\, m_A/M - \mathbf{p}_{\mathrm{rel}}$ and $\mathbf{p}_B = \mathbf{P}_{\mathrm{CoM}}\, m_B/M + \mathbf{p}_{\mathrm{rel}}$ into $d\mathbf{p}_A \wedge d\mathbf{R}_A + d\mathbf{p}_B \wedge d\mathbf{R}_B$, expand using $d\mathbf{R}_A = d\mathbf{R}_{\mathrm{CoM}} - (m_B/M)\, d\boldsymbol{\rho}$, $d\mathbf{R}_B = d\mathbf{R}_{\mathrm{CoM}} + (m_A/M)\, d\boldsymbol{\rho}$, and collect terms; the cross-terms cancel and the result is $d\mathbf{P}_{\mathrm{CoM}} \wedge d\mathbf{R}_{\mathrm{CoM}} + d\mathbf{p}_{\mathrm{rel}} \wedge d\boldsymbol{\rho}$. The 2-form is preserved.

**Lagrangian-side ($TQ$) corollary: kinetic metric block-diagonalises.** On the tangent bundle, the velocity transformation uses the *same* Jacobian $\mathbf{J}$ as the configuration-side transformation (linear point transformations have constant Jacobians, so the velocity and configuration sides transform identically). The kinetic-energy metric $\mathbf{g}_{\mathrm{old}} = \mathrm{diag}(m_A, m_B) \otimes \mathbb{I}_3$ transforms to:

$$
\mathbf{g}_{\mathrm{new}}
\;=\; (\mathbf{J}^{-1})^T\, \mathbf{g}_{\mathrm{old}}\, \mathbf{J}^{-1}
\;=\; \mathrm{diag}(M, \mu) \otimes \mathbb{I}_3.
$$

Direct verification (with the scalar block form of $\mathbf{J}^{-1}$ above): the $(1,2)$ and $(2,1)$ off-blocks vanish identically, and the diagonal blocks come out $M$ and $m_A m_B / M = \mu$ respectively. This is the **algebraic content** of the no-cross-term structure of $T_{\mathrm{CoM}}$ in §2.6: the CoM-decomposition Jacobian is the orthogonalising transformation of the kinetic metric with respect to the canonical $(m_A, m_B)$ inner product. It is *not* an accident; it is forced by the geometry of the centre-of-mass concept.

**Reduction perspective.** Because the Lagrangian $L = T - V$ depends on $\mathbf{R}_A, \mathbf{R}_B$ only through their difference $\boldsymbol{\rho}$ (translational invariance, Chapter 1 §1.6), the new coordinate $\mathbf{R}_{\mathrm{CoM}}$ is cyclic: $\partial L / \partial \mathbf{R}_{\mathrm{CoM}} = 0$. The conjugate momentum $\mathbf{P}_{\mathrm{CoM}}$ is therefore conserved (Noether's theorem applied to translation symmetry, as anticipated in Chapter 1 Table-of-Symmetries). The CoM coordinate can be eliminated entirely by **symplectic reduction by the translation subgroup** $\mathbb{R}^3 \subset SE(3)$ acting on $T^*(\mathbb{R}^3 \times \mathbb{R}^3)$, reducing the translational phase space from dimension 12 to dimension 6 (the relative motion). The general theory is in Marsden & Ratiu 1999 §1.2 + §11; for our purposes, fixing $\mathbf{P}_{\mathrm{CoM}} = \mathbf{0}$ and $\mathbf{R}_{\mathrm{CoM}} = \mathbf{0}$ (the CoM-rest-frame convention of Chapter 1 §1.3.1) accomplishes the reduction explicitly and leaves only the relative-motion variables $(\boldsymbol{\rho}, \mathbf{p}_{\mathrm{rel}})$ as active translational degrees of freedom in the remaining derivations.

> **Reference anchor.** Marsden & Ratiu 1999 §6.3 (cotangent lift of diffeomorphisms — the underlying theorem); §1.2 + §11 (symplectic reduction by group action — the framework that turns the CoM symmetry into the dimension count). Goldstein 2002 §9.4 covers canonical transformations in the classical-mechanics tradition without the modern geometric framing; useful as a coordinate-component cross-check. Arnold 1989 §3.16 is the concise modern treatment. The specific block-diagonalisation of the two-body kinetic metric is in Landau-Lifshitz Vol I §13 (reduced-mass framing, in coordinates).

## §2.7 Kinetic energy as a Riemannian metric on $TQ$ (geometric viewpoint)

Up to this point the chapter has worked in coordinates. The intrinsic / coordinate-free statement of what we have built is the following.

A **Riemannian metric** on a smooth manifold $M$ is a smooth assignment, to each tangent space $T_p M$, of a positive-definite symmetric bilinear form $g_p: T_p M \times T_p M \to \mathbb{R}$. The kinetic energy of a mechanical system whose configuration space is $M$ defines such a metric via $T = \tfrac{1}{2}\, g(\dot{q}, \dot{q})$. In coordinates, this becomes the familiar quadratic form $T = \tfrac{1}{2}\, g_{ij}(q)\, \dot{q}^i\, \dot{q}^j$.

For our system, $M = Q = SE(3)_A \times SE(3)_B$, and the metric $g$ is the one defined by the kinetic energy in §2.6. The reduced-mass / inertia structure of $g$ in the centre-of-mass frame is:

| Tangent subspace | Metric component |
|---|---|
| Relative-position velocity $\dot{\boldsymbol{\rho}}$ | $\mu\, \mathbb{I}_3$ (Euclidean, scaled by reduced mass) |
| Body $A$ angular velocity $\boldsymbol{\Omega}_A$ | $\mathbf{I}_A$ (body-frame inertia tensor) |
| Body $B$ angular velocity $\boldsymbol{\Omega}_B$ | $\mathbf{I}_B$ (body-frame inertia tensor) |
| Cross-blocks | All zero |

The metric is **block-diagonal** because the kinetic energy has no cross-terms (§2.6). The translational block is **flat** (Euclidean). The rotational blocks are not flat — they live on $SO(3)$, which is a Lie group with a non-trivial intrinsic geometry, and the metric on each $SO(3)$ factor is **left-invariant**: the metric value $g_q(\dot{q}, \dot{q})$ depends only on the body-frame angular velocity $\boldsymbol{\Omega}$, not on the current orientation $q$. This left-invariance is the geometric content of "the inertia tensor in the body frame is constant in time."

Left-invariant metrics on Lie groups are the source of structure that **Euler-Poincaré reduction** (Chapter 04) exploits. The Lagrangian for our system, when $T$ has this left-invariant form, reduces to equations of motion on the Lie algebra $\mathfrak{so}(3) \oplus \mathfrak{so}(3) \oplus \mathbb{R}^3$ (the relative-position tangent space) — a six-dimensional vector space — rather than on the 12-dimensional manifold $Q$ directly. This reduction is what produces **Euler's equations** for rigid-body rotation in the body frame as a structural consequence of the variational principle, rather than as a coordinate trick.

> **Reference anchor for the intrinsic perspective.** Arnold 1989 §3.27 (rigid body as a Lie group); Marsden & Ratiu 1999 §9.1 (left-invariant metrics on Lie groups); Marsden & Ratiu §13 (Euler-Poincaré reduction, anticipated here, derived in Chapter 04). For readers preferring a more elementary route, Goldstein §5.3 covers the inertia-tensor / kinetic-energy material in coordinates and Landau-Lifshitz Vol I §32 in concise form; both reach the same coordinate expressions as §2.6 without the geometric framing.

### §2.7.1 Left-invariance of the rotational kinetic-energy metric on $SO(3)$: explicit construction

The §2.7 main text asserted that the kinetic-energy metric on each $SO(3)$ factor of $Q$ is **left-invariant** and identified this as the geometric content of the body-frame inertia tensor's time-independence. The present subsection makes that statement concrete — what does it mean, mechanically, that "the metric is left-invariant" — and proves it via the explicit identification $g_R(\dot{q}_1, \dot{q}_2) = \langle \boldsymbol{\Omega}_1, \mathbf{I}\, \boldsymbol{\Omega}_2 \rangle_{\mathbb{R}^3}$ requested in OQ-2.3.

**Left translation on $SO(3)$ and its differential.** For each $h \in SO(3)$, **left translation** by $h$ is the diffeomorphism
$$
L_h : SO(3) \to SO(3), \qquad L_h(R) \;=\; h R.
$$
Its differential at $R$ is a linear map $(dL_h)_R : T_R SO(3) \to T_{hR} SO(3)$ between tangent spaces. For a tangent vector $\dot{q} \in T_R SO(3)$ — represented as the velocity of a curve through $R$ — the action of the differential is simply matrix multiplication on the left:
$$
(dL_h)_R\, \dot{q} \;=\; h \dot{q},
$$
which is again a tangent vector, now anchored at $h R$.

**Left-invariant vector fields.** A smooth vector field $X$ on $SO(3)$ is **left-invariant** if $(dL_h)_R X(R) = X(hR)$ for all $h, R \in SO(3)$. Equivalently, $X$ is determined throughout $SO(3)$ by its value $X(\mathbb{I}_3) \in T_{\mathbb{I}_3} SO(3) = \mathfrak{so}(3)$ at the identity, and propagated to other points by left translation. The space of left-invariant vector fields is therefore isomorphic to the **Lie algebra** $\mathfrak{so}(3)$ of antisymmetric $3 \times 3$ matrices, which in turn is isomorphic to $\mathbb{R}^3$ via the hat map of Chapter 1, §1.4. Both isomorphisms are linear and three-dimensional; together they justify treating angular-velocity vectors $\boldsymbol{\Omega} \in \mathbb{R}^3$ as the natural coordinates on $\mathfrak{so}(3)$.

**Left trivialisation of a tangent vector.** For any $\dot{q} \in T_R SO(3)$, the **left trivialisation** of $\dot{q}$ is the element $\widehat{\boldsymbol{\Omega}} \in \mathfrak{so}(3)$ obtained by pulling $\dot{q}$ back to the identity via $L_{R^{-1}}$:
$$
\widehat{\boldsymbol{\Omega}} \;\equiv\; (dL_{R^{-1}})_R\, \dot{q} \;=\; R^{-1} \dot{q} \;=\; R^T \dot{q}.
$$
Equivalently, $\dot{q} = R \widehat{\boldsymbol{\Omega}}$. This is exactly the body-frame relation from Chapter 1, §1.4: the body-frame angular velocity $\boldsymbol{\Omega}$ is the left-trivialisation of the orientation velocity $\dot{q}$. The body frame is the **left-trivialised frame** on $SO(3)$.

**The kinetic-energy metric is left-invariant.** Define the metric on $SO(3)$ by:
$$
\boxed{
g_R(\dot{q}_1, \dot{q}_2) \;\equiv\; \langle \boldsymbol{\Omega}_1, \mathbf{I}\, \boldsymbol{\Omega}_2 \rangle_{\mathbb{R}^3}, \qquad \widehat{\boldsymbol{\Omega}_i} = R^T \dot{q}_i,
}
$$
where $\mathbf{I}$ is the **constant** body-frame inertia tensor (a symmetric positive-definite $3 \times 3$ matrix, independent of $R$) and $\langle \cdot , \cdot \rangle_{\mathbb{R}^3}$ is the standard Euclidean inner product on $\mathbb{R}^3$.

To prove $g$ is left-invariant, evaluate at the translated point $h R$ with translated tangent vectors $(dL_h)_R \dot{q}_i = h \dot{q}_i$:
$$
g_{hR}\big((dL_h) \dot{q}_1,\, (dL_h) \dot{q}_2\big)
\;=\; \langle \boldsymbol{\Omega}'_1, \mathbf{I}\, \boldsymbol{\Omega}'_2 \rangle_{\mathbb{R}^3},
$$
where $\widehat{\boldsymbol{\Omega}'_i} = (hR)^T (h \dot{q}_i) = R^T h^T h \dot{q}_i = R^T \dot{q}_i = \widehat{\boldsymbol{\Omega}_i}$, using $h^T h = \mathbb{I}_3$ for $h \in SO(3)$. The left-trivialisations are unchanged by the left translation, hence:
$$
g_{hR}\big((dL_h) \dot{q}_1,\, (dL_h) \dot{q}_2\big)
\;=\; \langle \boldsymbol{\Omega}_1, \mathbf{I}\, \boldsymbol{\Omega}_2 \rangle_{\mathbb{R}^3}
\;=\; g_R(\dot{q}_1, \dot{q}_2). \qquad \blacksquare
$$

**Mechanical interpretation: the statement "$\mathbf{I}_{\mathrm{body}}$ is constant" *is* the statement "$g$ is left-invariant".** These are not two separate properties. The body-frame inertia tensor $\mathbf{I}$ being constant in time means it does not depend on the orientation $R$; equivalently, when the metric is expressed in the body-frame ($\boldsymbol{\Omega}$) parametrisation of tangent vectors, no $R$-dependent coefficients appear. The metric value $g_R(\dot{q}_1, \dot{q}_2) = \langle \boldsymbol{\Omega}_1, \mathbf{I}\, \boldsymbol{\Omega}_2 \rangle$ depends on $R$ only through the trivialisation $\boldsymbol{\Omega}_i = R^T \dot{q}_i$, and this dependence cancels out under simultaneous left-translation of $R$ and the tangent vectors. The cancellation is the algebraic content of left-invariance.

**Contrast: right trivialisation and inertial-frame description.** The mirror construction uses **right translation** $R_h(R) = R h$ instead of left. Its differential is right-multiplication: $(dR_h)_R \dot{q} = \dot{q} h$. The **right trivialisation** of $\dot{q}$ extracts $\widehat{\boldsymbol{\omega}} = \dot{q} R^{-1}$, the **inertial-frame** angular velocity from Chapter 1, §1.4. The same metric expressed in this trivialisation:
$$
g_R(\dot{q}_1, \dot{q}_2) \;=\; \langle \boldsymbol{\omega}_1, \mathbf{I}_R^{\mathrm{inertial}}\, \boldsymbol{\omega}_2 \rangle_{\mathbb{R}^3}, \qquad \mathbf{I}_R^{\mathrm{inertial}} = R\, \mathbf{I}\, R^T,
$$
uses the $R$-dependent **inertial-frame inertia tensor** $\mathbf{I}_R^{\mathrm{inertial}}$ — exactly the time-dependent tensor whose evolution we studied in §2.5.1. This same metric is now manifestly **right-invariant** (the right-translation analog of left-invariance, with the inertial frame playing the role the body frame played above), but it is *not* left-invariant: changing $R$ changes the matrix $\mathbf{I}_R^{\mathrm{inertial}}$. The same physical metric admits two complementary trivialisations; the *body* frame is the **left**-invariant frame and the *inertial* frame is the **right**-invariant frame.

This is the geometric source of an asymmetry already observed in §2.5: the inertia tensor is constant in the body frame (left-invariant frame) but the angular momentum is conserved in the inertial frame (right-equivariant). Constancy of $\mathbf{I}_{\mathrm{body}}$ comes from left-invariance of the metric; conservation of inertial-frame $\mathbf{L}$ comes from the right action (rotation of inertial frame) being a Noether symmetry of the system Lagrangian.

**Why left-invariance matters: preview of Euler-Poincaré reduction.** A Lagrangian that is left-invariant on a Lie group $G$ — i.e., $L(g, \dot{g}) = L(\mathbb{I}, g^{-1} \dot{g})$ for all $g$ — reduces to a Lagrangian $\ell(\boldsymbol{\Omega})$ on the Lie algebra $\mathfrak{g} = T_{\mathbb{I}} G$, where $\boldsymbol{\Omega} = g^{-1} \dot{g}$ is the left-trivialisation. The Euler-Lagrange equations on $TG$ then transform into the **Euler-Poincaré equations** on $\mathfrak{g}^*$:
$$
\frac{d}{dt} \frac{\partial \ell}{\partial \boldsymbol{\Omega}} \;=\; \mathrm{ad}^*_{\boldsymbol{\Omega}}\, \frac{\partial \ell}{\partial \boldsymbol{\Omega}}.
$$
For the rotational kinetic energy $\ell(\boldsymbol{\Omega}) = \tfrac{1}{2} \boldsymbol{\Omega}^T \mathbf{I}\, \boldsymbol{\Omega}$ on $\mathfrak{so}(3) \cong \mathbb{R}^3$, this specialises (with $\partial \ell / \partial \boldsymbol{\Omega} = \mathbf{I}\, \boldsymbol{\Omega}$ and $\mathrm{ad}^*_{\boldsymbol{\Omega}}\, v = -\boldsymbol{\Omega} \times v$ for $\mathfrak{so}(3)$ identified with $\mathbb{R}^3$) to:
$$
\frac{d}{dt}(\mathbf{I}\, \boldsymbol{\Omega}) \;=\; -\boldsymbol{\Omega} \times (\mathbf{I}\, \boldsymbol{\Omega}),
$$
which rearranges to **Euler's equations** $\mathbf{I}\, \dot{\boldsymbol{\Omega}} + \boldsymbol{\Omega} \times (\mathbf{I}\, \boldsymbol{\Omega}) = \mathbf{0}$ (torque-free; with external torque on the right-hand side in the general case). The body-frame Euler equations are therefore a **structural consequence** of the kinetic-energy metric's left-invariance, not a coordinate trick. Chapter 04 carries out this reduction in detail; Chapter 05 specialises to the body-frame Euler equations explicitly.

**Engineer-tier note.** The first-cut implementation `dynamics.py` parametrises orientation by a unit quaternion (a chart on $SO(3)$) and propagates the body-frame angular velocity $\boldsymbol{\Omega}_A$ via Euler's equations (lines 162–163), then reconstructs the orientation via the quaternion kinematics (lines 167–168). This is the **left-trivialised** integration: the dynamical variables $(\boldsymbol{\Omega}_A, q_A)$ are the left-trivialisation $\boldsymbol{\Omega}$ paired with the orientation $q$, and the equations of motion in these variables are time-autonomous (because $\mathbf{I}_A$ is constant). The alternative right-trivialised integration would propagate inertial-frame $\boldsymbol{\omega}_A$ with the time-dependent $\mathbf{I}_A^{\mathrm{inertial}}(t)$ — equivalent in physical content but algebraically more expensive per step, requiring the orientation-dependent inertia matrix to be updated at every RK4 sub-step. Left-trivialisation is the right choice for computational efficiency precisely because the metric is left-invariant.

> **Reference anchor.** Marsden & Ratiu 1999 §9.1 — the canonical exposition of left-invariant metrics on Lie groups, with $SO(3)$ as the worked example and Euler-Poincaré reduction in §13 as the application. Arnold 1989 §3.27 (rigid body as a Lie group) and Arnold appendix on Hamiltonian formalism cover the same material. For a tighter modern presentation, Holm-Schmah-Stoica 2009 *Geometric Mechanics and Symmetry* §6 carries the construction from $SO(3)$ to $SE(3)$ to general matrix Lie groups, useful preparation for Chapter 04's reduction of the full two-body Lagrangian.

> **Supplementary anchor for Czech-reading audiences.** Podolský 2018 Kapitola 6 §6.3 sets up the body-frame angular-velocity tensor via Eq. (6.4) $\Omega'_{ij} \equiv \Omega(\mathbf{e}'_i, \mathbf{e}'_j)$ — the components of the angular-velocity tensor evaluated on body-fixed (korotující) basis vectors. The trivialisation algebra at Eqs. (6.2) $\boldsymbol{\Omega}' = \frac{dA}{dt} A^t$ (body components) and (6.10) $\boldsymbol{\Omega} = A^t \frac{dA}{dt}$ (spatial components) is the engineer-tier rendering of the left-vs-right trivialisation distinction made geometrically here: Podolský's "$\boldsymbol{\Omega}'$ in korotující basis" is the left-trivialised body-frame angular velocity, his "$\boldsymbol{\Omega}$ in spatial basis" is the right-trivialised inertial-frame angular velocity, and the equivalence of the two physical descriptions corresponds to the same metric admitting both left- and right-invariant rendering on different sides of $SO(3)$. Podolský does not invoke Lie-group language explicitly, but the algebra he derives is the body/spatial trivialisation in coordinates and is what the modern Marsden-Ratiu treatment abstracts to coordinate-free form. Native-language anchor for the CS canonical-tier translation.

## §2.8 What is established by the end of this chapter

A reader who has worked through §2.1–§2.7 has:

- The kinetic energy $T$ of the two-rigid-body system as an explicit function on the tangent bundle $TQ$, derived from first principles by integration over the mass distribution of each body.
- König's decomposition into translational + rotational parts, with the centre-of-mass theorem as the algebraic basis for why no cross-term arises.
- The inertia tensor $\mathbf{I}_A$ as a symmetric positive-definite $3 \times 3$ matrix with the parallel-axis composition rule, mapping directly to the `geometries.py` composite-body construction in the first-cut engineer-tier code.
- The body-frame vs. inertial-frame angular velocity bookkeeping, with explicit identification of which frame keeps the inertia tensor static (body) vs. which frame the angular momentum is conserved in absence of external torque (inertial).
- The **rotational transport equation** $d\mathbf{I}^{\mathrm{inertial}}/dt = [\widehat{\boldsymbol{\omega}}, \mathbf{I}^{\mathrm{inertial}}]$ and its off-diagonal-element sensitivity to body symmetry and angular velocity (§2.5.1) — the engineer-tier bridge that translates the canonical formalism into practical sensitivity statements about initial conditions, body anisotropy, and integrator step-size selection.
- The reduced-mass form of translational kinetic energy in the centre-of-mass frame; the total $T_{\mathrm{CoM}}$ as the entry point to the Lagrangian.
- The **centre-of-mass decomposition as a symplectic point transformation** (§2.6.1) — the change of variables $(\mathbf{R}_A, \mathbf{R}_B) \to (\mathbf{R}_{\mathrm{CoM}}, \boldsymbol{\rho})$ given rigorous canonical-mechanics status as a cotangent-lift of a linear point transformation with unit-determinant Jacobian, with the induced momentum transformation $\mathbf{P}_{\mathrm{CoM}} = \mathbf{p}_A + \mathbf{p}_B$, $\mathbf{p}_{\mathrm{rel}} = \mu \dot{\boldsymbol{\rho}}$ and the kinetic-metric block-diagonalisation $\mathrm{diag}(M, \mu)$.
- The intrinsic geometric statement of $T$ as a Riemannian metric on $TQ$, block-diagonal, left-invariant on the $SO(3)$ factors. This is the structural property that Euler-Poincaré reduction (Chapter 04) will exploit.
- The **explicit construction** of left-invariance for the rotational kinetic-energy metric on each $SO(3)$ factor (§2.7.1), with the formula $g_R(\dot{q}_1, \dot{q}_2) = \langle \boldsymbol{\Omega}_1, \mathbf{I}\, \boldsymbol{\Omega}_2 \rangle$ derived and proven left-invariant by direct calculation. The mechanical equivalence "$\mathbf{I}_{\mathrm{body}}$ constant $\Leftrightarrow$ $g$ left-invariant" is identified; the contrast with the inertial-frame (right-trivialised) parametrisation is made explicit; Euler-Poincaré reduction is previewed as the consequence of left-invariance that produces Euler's equations in Chapter 04+05.

**What is NOT yet established:**

- The potential energy $V$ explicitly. This requires the mutual gravitational interaction between extended bodies, derived via multipole expansion. **Chapter 03.**
- The full Lagrangian $L = T - V$ and Euler-Lagrange equations of motion on $TQ$. **Chapter 04.**
- The body-frame Euler equations $\mathbf{I}\dot{\boldsymbol{\Omega}} + \boldsymbol{\Omega} \times (\mathbf{I}\boldsymbol{\Omega}) = \boldsymbol{\tau}$ as a specialisation of the Euler-Lagrange equations via Euler-Poincaré reduction. **Chapter 05.**
- Validation against analytical limits — torque-free symmetric top, asymmetric top instability, point-mass Kepler limit, etc. **Chapter 06**, with full catalogue locked in `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md`.

---

## §2.9 Open questions in this chapter

- **OQ-2.1.** ***Resolved 2026-05-23*** — addressed at new subsection §2.6.1 ("The centre-of-mass decomposition as a symplectic point transformation"). The change of variables $(\mathbf{R}_A, \mathbf{R}_B) \to (\mathbf{R}_{\mathrm{CoM}}, \boldsymbol{\rho})$ is the linear point transformation with Jacobian $\mathbf{J}$ ($\det \mathbf{J} = 1$); its cotangent lift to $T^*(\mathbb{R}^3 \times \mathbb{R}^3)$ is the symplectic transformation that maps $(\mathbf{p}_A, \mathbf{p}_B)$ to $(\mathbf{P}_{\mathrm{CoM}} = \mathbf{p}_A + \mathbf{p}_B,\, \mathbf{p}_{\mathrm{rel}} = \mu \dot{\boldsymbol{\rho}})$; the canonical 2-form is preserved (verified directly); the kinetic-energy metric block-diagonalises into $\mathrm{diag}(M, \mu) \otimes \mathbb{I}_3$. Formal anchor: Marsden & Ratiu 1999 §6.3 (cotangent lift) + §1.2 + §11 (symplectic reduction). On the tangent-bundle side, the velocity change of variables uses the same Jacobian $\mathbf{J}$ and produces the same block-diagonalisation of the kinetic metric.
- **OQ-2.2.** ***Resolved 2026-05-23*** — project-owner-supplied T&B PDF reviewed (1552-pp electronic edition). **Negative finding for König's theorem (§2.3):** the König-decomposition derivation does not appear in T&B's *Modern Classical Physics* volume; the rigid-body mechanics topic belongs to a companion mechanics volume not bundled with the uploaded edition. The §2.3 reference anchor remains Landau-Lifshitz Vol I §32. **Partial finding for inertia tensor (§2.4):** T&B §1.1.1 + §1.3 (coordinate-free geometric framing $\mathbf{I}(\,\cdot\,, \boldsymbol{\omega}) = \mathbf{J}$) added to the §2.4 reference block as the geometric-viewpoint complement to the coordinate derivation. The mass-distribution-integral definition of $\mathbf{I}$ in §2.4 retains its original anchors on Goldstein/LL/Hughes.
- **OQ-2.3.** ***Resolved 2026-05-23*** — addressed at new subsection §2.7.1 ("Left-invariance of the rotational kinetic-energy metric on $SO(3)$: explicit construction"). Concrete deliverable: left translation $L_h$ and its differential defined explicitly; left-invariant vector fields identified with $\mathfrak{so}(3) \cong \mathbb{R}^3$; left-trivialisation $\widehat{\boldsymbol{\Omega}} = R^T \dot{q}$ matched to the body-frame angular velocity convention of Chapter 1 §1.4; the metric $g_R(\dot{q}_1, \dot{q}_2) = \langle \boldsymbol{\Omega}_1, \mathbf{I}\, \boldsymbol{\Omega}_2 \rangle$ (boxed) shown left-invariant via direct algebraic computation; the mechanical equivalence "$\mathbf{I}_{\mathrm{body}}$ constant $\Leftrightarrow$ $g$ left-invariant" stated and proven; contrast with the right-trivialised / inertial-frame formulation given for context; preview of Euler-Poincaré reduction with the explicit Euler-equation specialisation; engineer-tier note tying the left-trivialised parametrisation directly to the `dynamics.py` Euler-equation integration choice. Anchors: Marsden & Ratiu 1999 §9.1 + §13; Arnold 1989 §3.27.
- **OQ-2.4.** ***Resolved 2026-05-23*** — engineer-tier bridge for the inertia-tensor formalism authored at new subsection §2.5.1 ("Time evolution of $\mathbf{I}_A^{\mathrm{inertial}}(t)$: off-diagonal sensitivity"). Captures the project-owner direction 2026-05-23: "Landau-Lifshitz is golden standard, yet for implementation of real calculations the time-derivative of the inertia tensor — especially its off-diagonal elements — would provide a bridge to illustrate how sensitive the dynamics are to initial mutual positions and how the axes and symmetries of the objects and the whole 2-body assembly enter." Deliverables: rotational transport equation $d\mathbf{I}^{\mathrm{inertial}}/dt = [\widehat{\boldsymbol{\omega}}, \mathbf{I}^{\mathrm{inertial}}]$ with explicit off-diagonal-rate formula $(d\mathbf{I}/dt)_{ij} = \varepsilon_{ijk}\, \omega_k\, (I_i - I_j)$ at $\mathbf{R}_q = \mathbb{I}$ instant; body-symmetry-class sensitivity table (spherical / axisymmetric / asymmetric); numerical sanity check for first-cut JWST-like body confirming $\tau_{\mathrm{off}} \approx 190\ \text{s}$ vs. integration window $600\ \text{s}$ and step size $0.05\ \text{s}$ → 4 orders of margin; explicit bridge to gravity-gradient-torque sensitivity for two-body assembly. Anchor: Hughes 2004 §4 + §5 + §8.

---

**Next chapter:** [03 — Mutual gravitational potential and multipole expansion](03-mutual-gravitational-potential.md) *(queued)*

**End 02-kinetic-energy.md**
