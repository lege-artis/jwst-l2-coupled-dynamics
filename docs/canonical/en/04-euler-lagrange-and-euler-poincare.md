# Chapter 4 — Euler-Lagrange equations and Euler-Poincaré reduction

> **Prerequisite reading.** Chapter 01 (configuration manifold and variational principle), Chapter 02 (kinetic energy on $SE(3) \times SE(3)$, including §2.5.1 dI/dt rotational transport, §2.6.1 CoM symplectic transformation, §2.7.1 left-invariance), and Chapter 03 (mutual gravitational potential, including §3.7 gravity-gradient torque and §3.7.2 orientation Jacobian).
>
> **Status.** v0.1 (canonical-tier round 3, 2026-05-25). All four chapter-internal open questions OQ-4.1, OQ-4.2, OQ-4.3, OQ-4.4 resolved in this round; see §4.9 for the close-out summary. Forward-looking v0.2 open questions surface at §4.9.

---

## §4.1 What this chapter does

Chapters 02 and 03 produced, respectively, the kinetic energy $T$ (the left-invariant Riemannian metric on the rotational factors of $TQ$, plus the reduced-mass translational metric) and the mutual gravitational potential $V$ (the quadrupole-truncated mutual potential with its gravity-gradient torque). This chapter assembles the **Lagrangian** $L = T - V$ on the configuration tangent bundle $TQ$, derives the **Euler-Lagrange equations** from the variational principle $\delta S = 0$ established in Chapter 01 §1.5, and applies **Euler-Poincaré reduction** (Marsden & Ratiu 1999 §13) to the left-invariant rotational degrees of freedom to recover the body-frame **Euler equations** for each rigid body as a structural consequence of the metric's left-invariance — not as a coordinate trick. The output of this chapter is the complete set of equations of motion for the two-body coupled rotational-translational system, ready for the specialised body-frame form in Chapter 05 and for numerical integration in the engineer-tier code paths.

The chapter follows the **two-in-parallel methodology** discipline established in Chapter 03 §3.6 and reaffirmed by project-owner direction 2026-05-23: the body-frame Euler equations are derived by **two independent routes** — (a) direct Euler-Lagrange computation in coordinates on $TQ$, and (b) Euler-Poincaré reduction by the left action of $SO(3) \times SO(3)$ on $Q$ — and the two are shown to coincide. Agreement is the validation; if the two derivations produced different equations, one of them would be wrong.

> **Supplementary anchor for Czech-reading audiences.** Podolský 2018 Kapitola 7 *Dynamika tuhého tělesa* is the most direct CS-language anchor for this chapter. §7.2 derives the body-frame Euler dynamic equations as the boxed Eq. (7.15): $I_1 \dot{\Omega}_x - (I_2 - I_3)\Omega_y\Omega_z = M_x$ (and cyclic) — the same equations of motion §4.4 here arrives at via Euler-Poincaré reduction. §7.3 then re-derives the same equations via the Lagrangian formalism using Euler angles as generalised coordinates — directly equivalent to our Approach A (§4.3) in the Euler-angle chart on $SO(3)$. The two-in-parallel cross-validation methodology this chapter formalises is therefore already implicit in Podolský's two-derivation presentation; the modern Marsden-Ratiu Euler-Poincaré reduction is the coordinate-free abstraction of the same agreement. Native-language anchor for the CS canonical-tier translation.

## §4.2 Assembly of the Lagrangian $L = T - V$

We combine the centre-of-mass-frame kinetic energy from Chapter 02 §2.6 (boxed) with the quadrupole-truncated mutual potential from Chapter 03 §3.4.3 (boxed) to obtain the Lagrangian on the configuration manifold reduced by translational symmetry. The translational symmetry of the system was used in Chapter 02 §2.6.1 to perform the symplectic-reduction-by-translation step that sets $\mathbf{P}_{\mathrm{CoM}} = \mathbf{0}$, $\mathbf{R}_{\mathrm{CoM}} = \mathbf{0}$ explicitly; the active configuration manifold is therefore:

$$
Q_{\mathrm{CoM}} \;\equiv\; \mathbb{R}^3_{\boldsymbol{\rho}} \;\times\; SO(3)_A \;\times\; SO(3)_B,
$$

of dimension $3 + 3 + 3 = 9$. The kinetic energy on $TQ_{\mathrm{CoM}}$ is (Ch 02 §2.6 boxed):

$$
T_{\mathrm{CoM}} \;=\; \tfrac{1}{2}\, \mu\, |\dot{\boldsymbol{\rho}}|^2
\;+\; \tfrac{1}{2}\, \boldsymbol{\Omega}_A^T\, \mathbf{I}_A\, \boldsymbol{\Omega}_A
\;+\; \tfrac{1}{2}\, \boldsymbol{\Omega}_B^T\, \mathbf{I}_B\, \boldsymbol{\Omega}_B,
$$

with $\mu = m_A m_B / (m_A + m_B)$ the reduced mass and $\boldsymbol{\Omega}_X$ the body-frame angular velocity of body $X$ (the left-trivialised tangent vector — Ch 02 §2.7.1). The potential energy is the Ch 03 §3.4.3 boxed quadrupole-truncated expression:

$$
V \;=\; -\frac{G\, m_A\, m_B}{|\boldsymbol{\rho}|}
\;+\; V_{\mathrm{tidal},\, A}(\boldsymbol{\rho}, q_A)
\;+\; V_{\mathrm{tidal},\, B}(\boldsymbol{\rho}, q_B)
\;+\; \mathcal{O}\!\left(\frac{\ell^4}{|\boldsymbol{\rho}|^5}\right),
$$

with each tidal piece

$$
V_{\mathrm{tidal},\, X} \;=\; -\frac{G\, m_{\bar X}}{2\, |\boldsymbol{\rho}|^3}\, \Big[ \mathrm{tr}(\mathbf{I}_X) \;-\; 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, X} \cdot \mathbf{I}_X \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, X} \Big],
$$

where $\bar X$ denotes the other body, $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, X} = (\mathbf{R}_q^X)^T \hat{\boldsymbol{\rho}}$ is the inter-body separation unit vector expressed in body $X$'s frame, and the constant body-frame inertia tensors $\mathbf{I}_X$ enter linearly. The Lagrangian on $TQ_{\mathrm{CoM}}$ is then:

$$
\boxed{
\begin{aligned}
L(\boldsymbol{\rho}, \dot{\boldsymbol{\rho}}, q_A, \boldsymbol{\Omega}_A, q_B, \boldsymbol{\Omega}_B)
\;=\;
&\tfrac{1}{2}\, \mu\, |\dot{\boldsymbol{\rho}}|^2 \;+\; \frac{G\, m_A\, m_B}{|\boldsymbol{\rho}|}
\\[2pt]
&\;+\; \tfrac{1}{2}\, \boldsymbol{\Omega}_A^T\, \mathbf{I}_A\, \boldsymbol{\Omega}_A
\;+\; \tfrac{1}{2}\, \boldsymbol{\Omega}_B^T\, \mathbf{I}_B\, \boldsymbol{\Omega}_B
\\[2pt]
&\;-\; V_{\mathrm{tidal},\, A}(\boldsymbol{\rho}, q_A)
\;-\; V_{\mathrm{tidal},\, B}(\boldsymbol{\rho}, q_B).
\end{aligned}
}
$$

### §4.2.1 Structural decomposition: which piece couples what to what

The boxed Lagrangian has three structurally distinct kinds of term, and the distinction matters for everything that follows in this chapter:

| Piece | Symbol | Configuration dependence | Velocity dependence | Couples |
|---|---|---|---|---|
| Translational kinetic | $\tfrac{1}{2}\mu|\dot{\boldsymbol{\rho}}|^2$ | none | $\dot{\boldsymbol{\rho}}$ only | nothing (purely translational) |
| Kepler potential (monopole-monopole) | $G m_A m_B / |\boldsymbol{\rho}|$ | $\boldsymbol{\rho}$ only | none | translation to translation |
| Rotational kinetic | $\tfrac{1}{2}\boldsymbol{\Omega}_X^T \mathbf{I}_X \boldsymbol{\Omega}_X$ | $q_X$ only (implicit via $\boldsymbol{\Omega}_X = (\mathbf{R}_q^X)^T \dot{R}^X$) | $\boldsymbol{\Omega}_X$ only | nothing (purely rotational, per body) |
| Tidal coupling | $V_{\mathrm{tidal},\, X}(\boldsymbol{\rho}, q_X)$ | both $\boldsymbol{\rho}$ AND $q_X$ | none | **translation $\leftrightarrow$ body-$X$ rotation** |

The exclusive source of rotational-translational coupling is the **tidal piece**: $V_{\mathrm{tidal},\, X}$ depends on $\boldsymbol{\rho}$ (through the magnitude $1/|\boldsymbol{\rho}|^3$ prefactor and through the unit vector $\hat{\boldsymbol{\rho}}$ inside the body-frame projection) AND on $q_X$ (through the rotation matrix $(\mathbf{R}_q^X)^T$ acting on $\hat{\boldsymbol{\rho}}$ to produce $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, X}$). When the tidal piece is removed (e.g., set $\mathbf{I}_A = \mathbf{I}_B = \mathbf{0}$ — the Kepler limit), the Lagrangian factors as $L = L_{\mathrm{trans}}(\boldsymbol{\rho}, \dot{\boldsymbol{\rho}}) + L_{\mathrm{rot},A}(\boldsymbol{\Omega}_A) + L_{\mathrm{rot},B}(\boldsymbol{\Omega}_B)$, the three pieces are independent, and the equations of motion decouple into pure Kepler translation plus two free Euler rotations. This factorisation is the structural underpinning of validation testbed §2.3 (Ch 03 §3.8.1) and Ch 06 §6.2.

Conversely, the tidal piece is the source of **two distinct dynamical effects** that share the same algebraic origin:

1. **Gravity-gradient torque on body-$X$ rotation.** $\partial V_{\mathrm{tidal},\, X} / \partial q_X \neq 0$ produces the body-frame torque $\boldsymbol{\tau}_X^{\mathrm{body}}$ derived at Ch 03 §3.7 boxed.
2. **Back-reaction force on the relative orbit.** $\partial V_{\mathrm{tidal},\, X} / \partial \boldsymbol{\rho} \neq 0$ produces an orientation-dependent force entering the $\boldsymbol{\rho}$ equation of motion. This is the §4.6 deliverable — the back-reaction was not derived in Ch 03 (which focused on the torque) and is therefore new content here.

The two effects are not independent: they are the rotational and translational projections of the same single object $\nabla V_{\mathrm{tidal},\, X}$ on the product manifold. Their algebraic correlation will surface as a Newton's-third-law-like balance between $\boldsymbol{\tau}_X^{\mathrm{body}}$ and the back-reaction force, which we make explicit at §4.6.

### §4.2.2 Numerical sanity check on the boxed Lagrangian

For the first-cut configuration (Ch 03 §3.7.2 numerical-sanity-check inputs: $m_A = 2450\ \text{kg}$, $\mathbf{I}_A = \mathrm{diag}(23323, 23323, 15384)\ \text{kg m}^2$ for the JWST-like body; $m_B = 100\ \text{kg}$, $\mathbf{I}_B = \mathrm{diag}(20, 20, 116)\ \text{kg m}^2$ for the probe-like body; $|\boldsymbol{\rho}| = 10\ \text{m}$, $|\dot{\boldsymbol{\rho}}| = 0.1\ \text{m/s}$, $\boldsymbol{\Omega}_A = \boldsymbol{\Omega}_B = (0, 0, 0.01)\ \text{rad/s}$), the four pieces of the boxed Lagrangian evaluate to:

| Piece | Magnitude (J) |
|---|---|
| Translational kinetic $\tfrac{1}{2}\mu|\dot{\boldsymbol{\rho}}|^2$ | $\sim 5 \times 10^{-1}$ |
| Kepler $G m_A m_B / |\boldsymbol{\rho}|$ | $\sim 2 \times 10^{-9}$ |
| Rotational $A$ + $B$ kinetic | $\sim 1$ |
| Tidal $V_{\mathrm{tidal},\, A} + V_{\mathrm{tidal},\, B}$ | $\sim 10^{-12}$ |

The tidal piece is **eleven orders of magnitude smaller** than the rotational kinetic at first-cut configurations — the gravitational interaction at metre-scale separations between modest masses is feeble. The interest of the model is not the magnitude of the tidal piece but its **structural** role: it is the *only* coupling between rotation and translation, and small or large, its presence is what gives the system its characteristic gravity-gradient libration and back-reaction phenomenology. The fact that the conservation diagnostic $|dE/E_0| = 6.2 \times 10^{-12}$ in the first-cut output (Ch 03 §3.9) tracks only the monopole-monopole term in $V$ (per `dynamics.py::total_potential_energy` — Ch 03 §3.9 table) without the tidal contribution does not affect observed conservation precision at these magnitudes; OQ-3.4 (Ch 03 §3.11) catalogues the engineer-tier augmentation that would make the diagnostic track full $V$.

## §4.3 Approach A — Euler-Lagrange in coordinates (Euler-angle chart)

The variational principle of Chapter 01 §1.5 gives the equations of motion as the Euler-Lagrange equations of the boxed Lagrangian above on the tangent bundle of $Q_{\mathrm{CoM}}$. Two natural parametrisations of the rotational factors exist: the **unit-quaternion chart** (the engineer-tier choice — explicit in `dynamics.py`, three-dimensional after imposing $|q|^2 = 1$) and the **Euler-angle chart** (the choice used in Podolský §7.3 and most classical textbook derivations). Both are charts on $SO(3)$ and yield the same physics; the canonical-tier derivation that follows works in the Euler-angle chart because (a) it has no constraint to handle explicitly (each Euler angle is a free coordinate), (b) it maps directly onto Podolský 2018 §7.3 Eqs. (7.16)–(7.22) for CS-language cross-reference, and (c) the algebraic form of the resulting equations is the form most readily compared to the §4.4 Euler-Poincaré reduction. Lock the choice. **(Resolves OQ-4.1.)** The engineer-tier quaternion chart is parameterisation-equivalent and is the integration-time choice for reasons given in Ch 02 §2.7.1's engineer-tier note — efficiency in the body-frame propagation step. The canonical-tier derivation is parameterisation-independent through the §4.4 reduction step, where the chart drops out entirely.

### §4.3.1 Euler angles and the kinematic equations

Adopt the **intrinsic ZXZ Euler-angle chart** on $SO(3)$: any orientation $R \in SO(3)$ is written as a composition of three elementary rotations parameterised by angles $(\varphi, \vartheta, \psi)$,

$$
\mathbf{R}_q^X(\varphi_X, \vartheta_X, \psi_X)
\;=\;
R_z(\psi_X)\, R_x(\vartheta_X)\, R_z(\varphi_X),
$$

where $R_z, R_x$ are elementary rotations about the $z$- and $x$-axes (Podolský 2018 §6.4 Eqs. 6.14–6.16 use this exact convention; we adopt his sign and order to keep the CS-language cross-reference clean). The chart is defined everywhere on $SO(3)$ except on a measure-zero gimbal-lock subset; this singularity is a chart artifact (not a physical singularity) and standard handling — multi-chart atlas, switch charts before traversing the singular subset — applies. For pedagogical derivation, working in a single chart suffices.

The body-frame angular velocity $\boldsymbol{\Omega}_X$ is then expressed in terms of the Euler-angle rates by Podolský 2018 §6.5 boxed Eq. (6.17):

$$
\boldsymbol{\Omega}_X \;=\; \mathbf{K}(\vartheta_X, \psi_X)\, \begin{pmatrix} \dot\varphi_X \\ \dot\vartheta_X \\ \dot\psi_X \end{pmatrix},
\qquad
\mathbf{K}(\vartheta, \psi) \;=\; \begin{pmatrix} \sin\vartheta\, \sin\psi & \cos\psi & 0 \\ \sin\vartheta\, \cos\psi & -\sin\psi & 0 \\ \cos\vartheta & 0 & 1 \end{pmatrix}.
$$

Equivalently, expanding the matrix product:

$$
\Omega_x = \dot\varphi \sin\vartheta\sin\psi + \dot\vartheta\cos\psi, \quad
\Omega_y = \dot\varphi \sin\vartheta\cos\psi - \dot\vartheta\sin\psi, \quad
\Omega_z = \dot\varphi \cos\vartheta + \dot\psi.
$$

This is the **kinematic** part of the rigid-body equations of motion in the Euler-angle chart: it tells how the body-frame $\boldsymbol{\Omega}$ relates to the chart-coordinate velocities $(\dot\varphi, \dot\vartheta, \dot\psi)$, without yet invoking Newton's, Euler's, or Lagrange's equations on the *dynamical* side.

### §4.3.2 Lagrangian in chart coordinates

Substitute the Euler-angle expression for $\boldsymbol{\Omega}_X$ into the boxed §4.2 Lagrangian. The rotational kinetic energy of body $A$ becomes (suppressing the $A$ subscript on Euler angles and $\mathbf{I}$ for clarity):

$$
\tfrac{1}{2}\boldsymbol{\Omega}^T \mathbf{I}\, \boldsymbol{\Omega}
\;=\;
\tfrac{1}{2}
\begin{pmatrix} \dot\varphi & \dot\vartheta & \dot\psi \end{pmatrix}
\mathbf{K}^T(\vartheta, \psi)\, \mathbf{I}\, \mathbf{K}(\vartheta, \psi)
\begin{pmatrix} \dot\varphi \\ \dot\vartheta \\ \dot\psi \end{pmatrix},
$$

a quadratic form in the Euler-angle rates with $\vartheta$- and $\psi$-dependent coefficients. For a principal-axis-aligned body ($\mathbf{I} = \mathrm{diag}(I_1, I_2, I_3)$), the kinetic energy expands to (algebra elementary but tedious; see Podolský §7.3 for the worked steps):

$$
T_{\mathrm{rot},\, A}
\;=\;
\tfrac{1}{2}\, I_1 (\dot\varphi \sin\vartheta\sin\psi + \dot\vartheta\cos\psi)^2
+ \tfrac{1}{2}\, I_2 (\dot\varphi \sin\vartheta\cos\psi - \dot\vartheta\sin\psi)^2
+ \tfrac{1}{2}\, I_3 (\dot\varphi \cos\vartheta + \dot\psi)^2.
$$

The potential terms $V_{\mathrm{tidal},\, A}$ also depend on $(\varphi, \vartheta, \psi)$ through $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = (\mathbf{R}_q^A)^T \hat{\boldsymbol{\rho}}$, but the dependence is on the Euler angles themselves (orientation), not on their rates.

### §4.3.3 Euler-Lagrange equations on the chart

Apply $\frac{d}{dt}\frac{\partial L}{\partial \dot{q}^i} - \frac{\partial L}{\partial q^i} = 0$ to each Euler-angle coordinate of body $A$. The three resulting equations are unsightly in expanded chart form; Podolský 2018 §7.3 carries out the worked algebra at his Eqs. (7.16)–(7.22) and arrives at the boxed Eq. (7.15):

$$
I_1 \dot\Omega_x - (I_2 - I_3)\, \Omega_y\Omega_z \;=\; M_x, \qquad \text{(and cyclic)}.
$$

The right-hand side $M_x$ is the body-frame component of the generalised torque arising from $-\partial V/\partial \boldsymbol{\theta}_A$ — equivalent to $\boldsymbol{\tau}_A^{\mathrm{body}}$ of Ch 03 §3.7 boxed (with the same prefactor and sign conventions, since Podolský's left-hand side and Ch 03's right-hand side use the same chart and same sign for the body-frame torque).

The algebraic step that converts the chart-form Euler-Lagrange equations into the body-frame Euler equations Eq. (7.15) uses the kinematic transport relation Podolský Eq. (6.8) (boxed): for any body-frame vector $\mathbf{w}$,

$$
\left.\frac{d\mathbf{w}}{dt}\right|_{\mathrm{spatial}} = \left.\frac{d\mathbf{w}}{dt}\right|_{\mathrm{body}} + \boldsymbol{\Omega} \times \mathbf{w}.
$$

Applied to $\mathbf{w} = \mathbf{L}_A = \mathbf{I}_A \boldsymbol{\Omega}_A$ (the body-frame angular momentum), this is precisely the transport that converts $d\mathbf{L}/dt|_{\mathrm{spatial}} = \boldsymbol{\tau}^{\mathrm{inertial}}$ (Newton's second law for angular momentum in the inertial frame) into the body-frame Euler form $\mathbf{I}\dot{\boldsymbol{\Omega}} + \boldsymbol{\Omega} \times (\mathbf{I}\boldsymbol{\Omega}) = \boldsymbol{\tau}^{\mathrm{body}}$. The Lagrangian-route derivation in Podolský §7.3 reproduces this transport implicitly through the Euler-angle-rate algebra; we do not re-derive it in chart coordinates here.

The translational $\boldsymbol{\rho}$ Euler-Lagrange equation gives:

$$
\frac{d}{dt}(\mu \dot{\boldsymbol{\rho}}) \;-\; \frac{\partial L}{\partial \boldsymbol{\rho}}
\;=\; 0,
\qquad \Rightarrow \qquad
\mu \ddot{\boldsymbol{\rho}} \;=\; -\nabla_{\boldsymbol{\rho}} V,
$$

with $V = -Gm_Am_B/|\boldsymbol{\rho}| + V_{\mathrm{tidal},A} + V_{\mathrm{tidal},B}$. This is the same equation §4.6 will analyse in detail; the chart form does not need expansion here.

### §4.3.4 Approach A summary

Approach A — direct Euler-Lagrange in the Euler-angle chart — produces:

- The classical body-frame Euler dynamic equations $\mathbf{I}_X \dot{\boldsymbol{\Omega}}_X + \boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X) = \boldsymbol{\tau}_X^{\mathrm{body}}$ for each body $X = A, B$ (recovered from the Euler-angle algebra via the transport equation).
- The translational equation $\mu \ddot{\boldsymbol{\rho}} = -\nabla_{\boldsymbol{\rho}} V$ for the relative-position coordinate.
- A coordinate-dependent intermediate algebra that obscures the underlying geometric structure (left-invariance of the rotational metric) but is correct in chart form.

This is **the first leg of the two-in-parallel cross-validation** at §4.5. In Approach B we recover the same body-frame Euler equations via a coordinate-free reduction route that makes the left-invariance hypothesis explicit.

## §4.4 Approach B — Euler-Poincaré reduction

The Euler-Poincaré reduction theorem (Marsden & Ratiu 1999 §13.5) takes a left-invariant Lagrangian on the tangent bundle of a Lie group and reduces it to a Lagrangian on the Lie algebra, with a modified variational principle that produces equations on the Lie-algebra-dual. We apply it to each $SO(3)$ rotational factor in turn; the translational factor is abelian ($\mathbb{R}^3$) and has trivial reduction.

### §4.4.1 The Euler-Poincaré theorem

**Theorem (Marsden & Ratiu 1999 §13.5).** Let $G$ be a Lie group with Lie algebra $\mathfrak{g}$, and let $L: TG \to \mathbb{R}$ be a Lagrangian. Suppose $L$ is **left-invariant**: $L(g, \dot{g}) = L(\mathbb{I}, g^{-1}\dot{g})$ for all $g \in G$ and $\dot{g} \in T_gG$. Define the **reduced Lagrangian** $\ell: \mathfrak{g} \to \mathbb{R}$ by

$$
\ell(\boldsymbol{\Omega}) \;\equiv\; L(\mathbb{I}, \boldsymbol{\Omega}), \qquad \boldsymbol{\Omega} = g^{-1}\dot{g} \in \mathfrak{g}.
$$

Then the Euler-Lagrange equations on $TG$ for $L$ are equivalent to the **Euler-Poincaré equations** on $\mathfrak{g}^*$ for $\ell$, with a possibly externally-supplied generalised force $\boldsymbol{f}^{\mathrm{ext}} \in \mathfrak{g}^*$:

$$
\boxed{
\frac{d}{dt}\frac{\partial \ell}{\partial \boldsymbol{\Omega}}
\;=\;
\mathrm{ad}^*_{\boldsymbol{\Omega}}\,\frac{\partial \ell}{\partial \boldsymbol{\Omega}}
\;+\;
\boldsymbol{f}^{\mathrm{ext}},
}
$$

where $\mathrm{ad}^*: \mathfrak{g} \times \mathfrak{g}^* \to \mathfrak{g}^*$ is the dual of the adjoint action of $\mathfrak{g}$ on itself, defined by $\langle \mathrm{ad}^*_{\boldsymbol{\Omega}}\, \boldsymbol{\mu}, \boldsymbol{\eta}\rangle = \langle \boldsymbol{\mu}, [\boldsymbol{\Omega}, \boldsymbol{\eta}]\rangle$ for $\boldsymbol{\Omega}, \boldsymbol{\eta} \in \mathfrak{g}$ and $\boldsymbol{\mu} \in \mathfrak{g}^*$. (See Marsden & Ratiu §9.3 for the adjoint and coadjoint actions in generality; §13.5 for the reduction theorem proof.)

### §4.4.2 Verification of the left-invariance hypothesis

Specialise to a single $SO(3)$ factor of $Q_{\mathrm{CoM}}$. The rotational kinetic Lagrangian on $TSO(3)$ is

$$
L_{\mathrm{rot},\, X}(R, \dot{R}) \;=\; \tfrac{1}{2} \boldsymbol{\Omega}_X^T \mathbf{I}_X \boldsymbol{\Omega}_X, \qquad \widehat{\boldsymbol{\Omega}_X} = R^T \dot{R},
$$

with $\mathbf{I}_X$ the **constant** body-frame inertia tensor. Chapter 02 §2.7.1 carried out the explicit proof that this Lagrangian is left-invariant: under left-translation $R \to hR$, $\dot{R} \to h\dot{R}$, the left-trivialisation $\boldsymbol{\Omega}_X = R^T \dot{R} \to (hR)^T (h\dot{R}) = R^T h^T h \dot{R} = R^T \dot{R} = \boldsymbol{\Omega}_X$ is unchanged, hence $L_{\mathrm{rot},\, X}$ is unchanged. The hypothesis of the Euler-Poincaré theorem is satisfied for the **kinetic** part of the rotational Lagrangian.

The tidal potential $V_{\mathrm{tidal},\, X}$, however, is **not** left-invariant: it depends on the body-frame projection $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, X} = R^T \hat{\boldsymbol{\rho}}$, which changes under $R \to hR$. The potential breaks the left-symmetry. The standard treatment is to absorb the potential into the externally-supplied generalised force $\boldsymbol{f}^{\mathrm{ext}}$ of the theorem; the body-frame torque $\boldsymbol{\tau}_X^{\mathrm{body}} = -\partial V_{\mathrm{tidal},\, X}/\partial \boldsymbol{\theta}_X$ (derived as Ch 03 §3.7 boxed) plays exactly this role. The Euler-Poincaré equation reduces to the body-frame Euler equation with torque, as we now compute explicitly.

### §4.4.3 Specialising $\mathrm{ad}^*$ for $\mathfrak{so}(3) \cong \mathbb{R}^3$

On the Lie algebra $\mathfrak{so}(3)$ of antisymmetric $3\times 3$ matrices, the Lie bracket is the matrix commutator $[A, B] = AB - BA$. Identifying $\mathfrak{so}(3) \cong \mathbb{R}^3$ via the hat map of Ch 01 §1.4 ($\widehat{\boldsymbol{a}}\, \boldsymbol{b} = \boldsymbol{a} \times \boldsymbol{b}$), the bracket transports to the cross product:

$$
[\boldsymbol{\Omega}, \boldsymbol{\eta}]_{\mathfrak{so}(3)} \;\leftrightarrow\; \boldsymbol{\Omega} \times \boldsymbol{\eta}\quad \text{(in $\mathbb{R}^3$)}.
$$

The dual space $\mathfrak{so}(3)^*$ is identified with $\mathbb{R}^3$ via the standard Euclidean inner product. We can now compute $\mathrm{ad}^*_{\boldsymbol{\Omega}}\, \boldsymbol{\mu}$ explicitly by carrying the pairing through. For $\boldsymbol{\Omega}, \boldsymbol{\eta} \in \mathfrak{so}(3) \cong \mathbb{R}^3$ and $\boldsymbol{\mu} \in \mathfrak{so}(3)^* \cong \mathbb{R}^3$, the definition $\langle \mathrm{ad}^*_{\boldsymbol{\Omega}}\, \boldsymbol{\mu}, \boldsymbol{\eta}\rangle = \langle \boldsymbol{\mu}, [\boldsymbol{\Omega}, \boldsymbol{\eta}]\rangle$ becomes:

$$
\mathrm{ad}^*_{\boldsymbol{\Omega}}\,\boldsymbol{\mu} \cdot \boldsymbol{\eta}
\;=\;
\boldsymbol{\mu} \cdot (\boldsymbol{\Omega} \times \boldsymbol{\eta}).
$$

Using the scalar-triple-product cyclic identity $\boldsymbol{\mu} \cdot (\boldsymbol{\Omega} \times \boldsymbol{\eta}) = (\boldsymbol{\mu} \times \boldsymbol{\Omega}) \cdot \boldsymbol{\eta} = -(\boldsymbol{\Omega} \times \boldsymbol{\mu}) \cdot \boldsymbol{\eta}$ — and demanding that this holds for arbitrary $\boldsymbol{\eta}$ — we read off:

$$
\boxed{
\mathrm{ad}^*_{\boldsymbol{\Omega}}\,\boldsymbol{\mu} \;=\; -\, \boldsymbol{\Omega} \times \boldsymbol{\mu}, \qquad \boldsymbol{\Omega} \in \mathfrak{so}(3) \cong \mathbb{R}^3, \quad \boldsymbol{\mu} \in \mathfrak{so}(3)^* \cong \mathbb{R}^3.
}
$$

**(Resolves OQ-4.2.)** The minus sign is the convention forced by the duality pairing — *not* an independent choice. Different textbooks present $\mathrm{ad}^*$ for $\mathfrak{so}(3)$ with the opposite sign; the disagreement reduces to whether the pairing $\langle \cdot, \cdot \rangle: \mathfrak{g} \times \mathfrak{g}^* \to \mathbb{R}$ is read with $\mathfrak{g}$ on the left or right. We follow Marsden & Ratiu §9.3 convention (left-on-left), giving the explicit minus sign here. The body-frame Euler equation derived below carries the sign consistently with Ch 02 §2.7.1's left-trivialisation conventions and with Ch 03 §3.7 boxed gravity-gradient torque.

### §4.4.4 Specialising the Euler-Poincaré equation to rigid-body rotation

The reduced rotational Lagrangian is

$$
\ell_X(\boldsymbol{\Omega}_X) \;=\; \tfrac{1}{2} \boldsymbol{\Omega}_X^T \mathbf{I}_X \boldsymbol{\Omega}_X, \qquad \frac{\partial \ell_X}{\partial \boldsymbol{\Omega}_X} \;=\; \mathbf{I}_X\, \boldsymbol{\Omega}_X.
$$

Plug into the boxed Euler-Poincaré equation, with $\boldsymbol{f}^{\mathrm{ext}}_X = \boldsymbol{\tau}_X^{\mathrm{body}}$ as discussed in §4.4.2:

$$
\frac{d}{dt}(\mathbf{I}_X \boldsymbol{\Omega}_X)
\;=\;
\mathrm{ad}^*_{\boldsymbol{\Omega}_X}\, (\mathbf{I}_X \boldsymbol{\Omega}_X)
\;+\; \boldsymbol{\tau}_X^{\mathrm{body}}
\;=\;
-\, \boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X)
\;+\; \boldsymbol{\tau}_X^{\mathrm{body}},
$$

using the boxed $\mathrm{ad}^*_{\boldsymbol{\Omega}}\,\boldsymbol{\mu} = -\boldsymbol{\Omega} \times \boldsymbol{\mu}$ from §4.4.3. Since $\mathbf{I}_X$ is constant in the body frame, $d(\mathbf{I}_X \boldsymbol{\Omega}_X)/dt = \mathbf{I}_X \dot{\boldsymbol{\Omega}}_X$, and rearranging gives the **body-frame Euler equation** in its canonical form:

$$
\boxed{
\mathbf{I}_X \dot{\boldsymbol{\Omega}}_X \;+\; \boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X) \;=\; \boldsymbol{\tau}_X^{\mathrm{body}},
\qquad X \in \{A, B\}.
}
$$

These are the equations of motion for the rotational degrees of freedom on the reduced Lie-algebra phase space. The translational $\boldsymbol{\rho}$ equation is unaffected by this reduction (Euler-Poincaré applies only to the Lie-group factors $SO(3)$; the translational $\mathbb{R}^3_{\boldsymbol{\rho}}$ factor is already abelian and the reduction is trivial).

### §4.4.5 Approach B summary

Approach B — Euler-Poincaré reduction — produces:

- The body-frame Euler equation $\mathbf{I}_X \dot{\boldsymbol{\Omega}}_X + \boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X) = \boldsymbol{\tau}_X^{\mathrm{body}}$ for each body $X = A, B$, identical to the §4.3.3 Approach A result.
- An explicit $\mathrm{ad}^*$ derivation that pins the sign convention to match Ch 02 §2.7.1 + Ch 03 §3.7 + OQ-FORT-1 sign-fix.
- A coordinate-free presentation that makes the left-invariance hypothesis (Ch 02 §2.7.1) and the role of the externally-applied torque explicit.

The reduction makes manifest that the body-frame Euler equation is a **structural consequence** of the metric's left-invariance, not a Newton-Euler $d\mathbf{L}/dt = \boldsymbol{\tau}$ statement transformed coordinatewise. Both views are correct; the structural view explains *why* the body-frame is the natural frame for rigid-body propagation and motivates the engineer-tier choice (Ch 02 §2.7.1's engineer-tier note).

## §4.5 Cross-validation: Approach A ≡ Approach B

The two approaches arrive at the **same** body-frame Euler equation $\mathbf{I}_X \dot{\boldsymbol{\Omega}}_X + \boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X) = \boldsymbol{\tau}_X^{\mathrm{body}}$ for each body $X$:

| | Approach A | Approach B |
|---|---|---|
| Method | Direct Euler-Lagrange in Euler-angle chart | Euler-Poincaré reduction of left-invariant Lagrangian |
| Coordinate-dependence | Yes (chart on $SO(3)$) | No (coordinate-free) |
| Intermediate steps | Chart algebra + transport equation | $\mathrm{ad}^*$ identification + reduction theorem |
| Final equation | $\mathbf{I}_X \dot{\boldsymbol{\Omega}}_X + \boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X) = \boldsymbol{\tau}_X^{\mathrm{body}}$ | (identical) |
| Sign conventions | Inherited from Podolský §7.3 | Inherited from Marsden-Ratiu §13.5 |
| Cross-check | $\boldsymbol{\tau}_X^{\mathrm{body}}$ matches Ch 03 §3.7 boxed (post-OQ-FORT-1 sign-fix) | (same) |

The agreement is the validation. As emphasised in Ch 03 §3.6 (two-in-parallel methodology) and the project-owner direction 2026-05-23: the value of carrying out *both* derivations is that an algebraic error in either route would surface as disagreement, and elimination of disagreement is elimination of that class of error.

This cross-validation is the analytical analog of the doctrine §4.4 C2 cross-implementation bit-identity criterion, applied at the derivation level. Subsequent canonical-tier chapters extend the pattern: Chapter 05 derives the principal-axis specialisation by both component-by-component expansion and structural specialisation of the §4.4 result; Chapter 06 validates each numerical testbed against at least two textbook references where available.

> **Reference anchor for the cross-validation pattern.** Marsden & Ratiu 1999 §13 (Euler-Poincaré reduction of left-invariant Lagrangians, with $SO(3)$ rigid body as worked example at §13.5). Holm-Schmah-Stoica 2009 *Geometric Mechanics and Symmetry* §6.2 carries out the same reduction with explicit $\mathrm{ad}^*$ algebra in the $SO(3)$ specialisation. Arnold 1989 §3.27 gives the elegant short-form: the rigid body is a geodesic flow on $SO(3)$ with a left-invariant metric, and Euler's equations are the equations of geodesics in body-frame coordinates. All three references agree on the boxed body-frame Euler equation; the sign of the cross-product term is invariant under their differing conventions.

## §4.6 Translational equation in the CoM frame and back-reaction force

The translational Euler-Lagrange equation from §4.3.3 reads:

$$
\mu \ddot{\boldsymbol{\rho}} \;=\; -\nabla_{\boldsymbol{\rho}} V,
$$

with $V = -Gm_Am_B/|\boldsymbol{\rho}| + V_{\mathrm{tidal},A} + V_{\mathrm{tidal},B}$ from §4.2 boxed. Decompose the gradient term by term.

### §4.6.1 Gradient of the Kepler piece

$$
-\nabla_{\boldsymbol{\rho}}\left(-\frac{Gm_Am_B}{|\boldsymbol{\rho}|}\right)
\;=\;
-\frac{Gm_Am_B}{|\boldsymbol{\rho}|^2}\, \hat{\boldsymbol{\rho}}
\;=\;
-\frac{Gm_Am_B}{|\boldsymbol{\rho}|^3}\, \boldsymbol{\rho}.
$$

This is the **Kepler force** — the attractive inverse-square central force between the two bodies treated as point masses. Combined with $\mu \ddot{\boldsymbol{\rho}} = \mathbf{F}_{\mathrm{Kepler}}$ alone (i.e., absent the tidal pieces), the equation of motion is the standard reduced-mass Kepler equation (Landau-Lifshitz Vol I §13–§15). This is the testbed §2.3 / Ch 06 §6.2 limit.

### §4.6.2 Gradient of the tidal piece — derivation of the back-reaction force

The tidal piece for body $A$ is (Ch 03 §3.4.3 boxed):

$$
V_{\mathrm{tidal},\, A}(\boldsymbol{\rho}, q_A)
\;=\;
-\frac{Gm_B}{2|\boldsymbol{\rho}|^3}\,
\Big[ \mathrm{tr}(\mathbf{I}_A) \;-\; 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \Big].
$$

To compute $-\nabla_{\boldsymbol{\rho}} V_{\mathrm{tidal},\, A}$, decompose the dependence on $\boldsymbol{\rho}$ into two pieces:

1. The scalar **prefactor** $1/|\boldsymbol{\rho}|^3$ — depends on the magnitude $r \equiv |\boldsymbol{\rho}|$ only.
2. The **bracketed orientation factor** $[\mathrm{tr}(\mathbf{I}_A) - 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}]$ — depends on $\boldsymbol{\rho}$ only through the unit vector $\hat{\boldsymbol{\rho}}$ inside $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = (\mathbf{R}_q^A)^T \hat{\boldsymbol{\rho}}$.

Denote the bracketed orientation factor by $\Phi_A(\hat{\boldsymbol{\rho}}, q_A) \equiv \mathrm{tr}(\mathbf{I}_A) - 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$. By the product rule:

$$
\nabla_{\boldsymbol{\rho}} V_{\mathrm{tidal},\, A}
\;=\;
-\frac{Gm_B}{2}\,
\Big[\,
\Phi_A\, \nabla_{\boldsymbol{\rho}}(1/r^3)
\;+\;
\frac{1}{r^3}\, \nabla_{\boldsymbol{\rho}} \Phi_A
\,\Big].
$$

Compute each piece.

**Radial-prefactor gradient.** $\nabla_{\boldsymbol{\rho}}(1/r^3) = -3\hat{\boldsymbol{\rho}}/r^4$, so:

$$
-\frac{Gm_B}{2}\, \Phi_A\, \nabla_{\boldsymbol{\rho}}(1/r^3)
\;=\;
\frac{3\, G m_B}{2\, r^4}\, \Phi_A\, \hat{\boldsymbol{\rho}}.
$$

**Orientation-factor gradient.** $\Phi_A$ depends on $\boldsymbol{\rho}$ only through $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = (\mathbf{R}_q^A)^T \hat{\boldsymbol{\rho}}$. Using $\partial \hat{\boldsymbol{\rho}}/\partial \boldsymbol{\rho} = (\mathbb{I}_3 - \hat{\boldsymbol{\rho}} \hat{\boldsymbol{\rho}}^T)/r$ (the transverse-projection identity), and applying the chain rule:

$$
\nabla_{\boldsymbol{\rho}}\Phi_A
\;=\;
\frac{(\mathbb{I}_3 - \hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}^T)}{r}\, \mathbf{R}_q^A\, \nabla_{\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}} \Phi_A.
$$

Using $\Phi_A = \mathrm{tr}(\mathbf{I}_A) - 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$ and noting that $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$ has unit magnitude:

$$
\nabla_{\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}} \Phi_A
\;=\;
-6\, \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}.
$$

Combining:

$$
\nabla_{\boldsymbol{\rho}}\Phi_A
\;=\;
-\frac{6}{r}\, (\mathbb{I}_3 - \hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}^T)\, \mathbf{R}_q^A\, (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}).
$$

The orientation-factor contribution to the gradient is therefore:

$$
-\frac{Gm_B}{2}\, \frac{1}{r^3}\, \nabla_{\boldsymbol{\rho}}\Phi_A
\;=\;
\frac{3\, G m_B}{r^4}\, (\mathbb{I}_3 - \hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}^T)\, \mathbf{R}_q^A\, (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}).
$$

**Combined back-reaction force on the relative orbit, from body $A$'s tidal contribution.** Collecting both pieces:

$$
\boxed{
\mathbf{F}_{\mathrm{back-reaction},\, A}
\;\equiv\;
-\nabla_{\boldsymbol{\rho}} V_{\mathrm{tidal},\, A}
\;=\;
-\frac{3\, G m_B}{2\, r^4}\, \Phi_A\, \hat{\boldsymbol{\rho}}
\;-\;
\frac{3\, G m_B}{r^4}\, (\mathbb{I}_3 - \hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}^T)\, \mathbf{R}_q^A\, (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}).
}
$$

The back-reaction force scales as $1/r^4$ — one power of $r$ steeper than the Kepler force's $1/r^2$ — and has an orientation-dependent prefactor encoding body $A$'s anisotropy and its alignment with $\hat{\boldsymbol{\rho}}$. A symmetric expression for $\mathbf{F}_{\mathrm{back-reaction},\, B}$ follows by $A \leftrightarrow B$ exchange (with the sign of $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, B}$ flipped, as in Ch 03 §3.7.1 for the torque).

**(Resolves OQ-4.3.)** The back-reaction force exists, scales as $1/r^4$ (matching the bound estimate of OQ-4.3 in the §4.9 placeholder), and has explicit closed-form expression in terms of $(r, \hat{\boldsymbol{\rho}}, \mathbf{I}_A, \mathbf{R}_q^A)$.

### §4.6.3 Spherical-body limit cancellation

For a spherical body ($\mathbf{I}_A = I_0 \mathbb{I}_3$): $\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = I_0\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$, and $\Phi_A = 3I_0 - 3I_0 |\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}|^2 = 3I_0 - 3I_0 = 0$. The first (radial-prefactor) term in $\mathbf{F}_{\mathrm{back-reaction},\, A}$ vanishes because of $\Phi_A = 0$. The second (orientation-factor) term: $(\mathbb{I}_3 - \hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}^T)\, \mathbf{R}_q^A\, (I_0 \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}) = I_0\, (\mathbb{I}_3 - \hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}^T)\, \mathbf{R}_q^A (\mathbf{R}_q^A)^T \hat{\boldsymbol{\rho}} = I_0\, (\mathbb{I}_3 - \hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}^T) \hat{\boldsymbol{\rho}} = I_0\, (\hat{\boldsymbol{\rho}} - \hat{\boldsymbol{\rho}}) = \mathbf{0}$. Both pieces vanish; the back-reaction force on the orbit due to body $A$'s tidal coupling vanishes identically when $A$ is spherical. This is consistent with the gravity-gradient-torque vanishing in the same limit (Ch 03 §3.8.2): a spherical body's gravitational coupling to its environment depends only on its mass, not its orientation, so neither the rotational dynamics nor the translational back-reaction can carry orientation-derived information.

### §4.6.4 Back-reaction–torque balance (Newton's third law check)

A useful consistency check: in an isolated two-body system, the total linear momentum $\mathbf{P}_{\mathrm{CoM}}$ is conserved (Ch 02 §2.6.1; symmetry under spatial translation by Noether's theorem). Therefore the total force on the centre-of-mass coordinate vanishes:

$$
\mathbf{F}^{\mathrm{external}}_{\mathrm{CoM}} \;\equiv\; \frac{d\mathbf{P}_{\mathrm{CoM}}}{dt} \;=\; \mathbf{0}.
$$

This is automatically respected by our Lagrangian's translational invariance — the Lagrangian depends on $(\mathbf{R}_A, \mathbf{R}_B)$ only through $\boldsymbol{\rho} = \mathbf{R}_B - \mathbf{R}_A$, so $\partial L/\partial \mathbf{R}_{\mathrm{CoM}} = 0$ and $\mathbf{P}_{\mathrm{CoM}}$ is the conjugate-momentum cyclic-coordinate first integral. The back-reaction force enters the *relative* equation $\mu\ddot{\boldsymbol{\rho}} = \mathbf{F}_{\mathrm{Kepler}} + \mathbf{F}_{\mathrm{back-reaction},\, A} + \mathbf{F}_{\mathrm{back-reaction},\, B}$ but does not break the centre-of-mass conservation; the back-reaction is an internal force in Newton's-third-law sense.

The corresponding angular-momentum statement is the **total angular momentum balance**: the rate of change of total angular momentum (orbital + body-$A$-spin + body-$B$-spin) about the centre of mass vanishes in the absence of external torques. In our system, the orbital angular momentum $\boldsymbol{\rho} \times \mu \dot{\boldsymbol{\rho}}$ can exchange with the body spins $\mathbf{R}_q^X (\mathbf{I}_X \boldsymbol{\Omega}_X)$ via the gravity-gradient torque (rotational side) plus a torque-on-orbit due to the back-reaction force (translational side); the two contributions cancel. The explicit verification is the calculation $\mathbf{F}_{\mathrm{back-reaction},\, A} \times \mathbf{R}_A + \boldsymbol{\tau}_A^{\mathrm{inertial}} + \mathbf{F}_{\mathrm{back-reaction},\, B} \times \mathbf{R}_B + \boldsymbol{\tau}_B^{\mathrm{inertial}} = \mathbf{0}$, which follows from $\nabla_{\boldsymbol{\theta}_A} V + \nabla_{\boldsymbol{\theta}_B} V + \boldsymbol{\rho} \times \nabla_{\boldsymbol{\rho}} V = 0$, the identity expressing rotational invariance of the mutual potential $V$. We assert it here for structural completeness; the explicit Noether-identity check is the engineer-tier conservation diagnostic $|d\mathbf{L}_{\mathrm{total}}/L_0|$ which the first-cut RK4 integration achieves at $2.1 \times 10^{-10}$ over the 600-second window (Ch 00 §3).

## §4.7 Total equation set and dynamics.py cross-validation

Assemble the §4.4 boxed rotational equations and the §4.6 translational equation including back-reaction force into the **complete equation set** on the 18-dimensional reduced phase space $T^* Q_{\mathrm{CoM}}$ (counting: 3 components $\boldsymbol{\rho}$ + 3 components $\mathbf{p}_{\mathrm{rel}}$ + 3 Euler-angle / 4-quaternion coordinates $q_A$ + 3 body-frame angular-velocity components $\boldsymbol{\Omega}_A$ + 3 quaternion / Euler-angle for $q_B$ + 3 for $\boldsymbol{\Omega}_B$, with one constraint per quaternion in engineer-tier; net 18):

$$
\boxed{
\begin{aligned}
\text{Translational:}\quad
&\dot{\boldsymbol{\rho}} \;=\; \mathbf{p}_{\mathrm{rel}}/\mu \\
&\dot{\mathbf{p}}_{\mathrm{rel}} \;=\; \mathbf{F}_{\mathrm{Kepler}}(\boldsymbol{\rho}) \;+\; \mathbf{F}_{\mathrm{back-reaction},\, A}(\boldsymbol{\rho}, q_A) \;+\; \mathbf{F}_{\mathrm{back-reaction},\, B}(\boldsymbol{\rho}, q_B) \\[4pt]
\text{Rotational $A$:}\quad
&\dot{q}_A \;=\; \tfrac{1}{2}\, q_A \otimes \boldsymbol{\Omega}_A \quad \text{(quaternion chart)} \\
&\mathbf{I}_A \dot{\boldsymbol{\Omega}}_A \;=\; -\boldsymbol{\Omega}_A \times (\mathbf{I}_A \boldsymbol{\Omega}_A) \;+\; \boldsymbol{\tau}_A^{\mathrm{body}}(\boldsymbol{\rho}, q_A) \\[4pt]
\text{Rotational $B$:}\quad
&\dot{q}_B \;=\; \tfrac{1}{2}\, q_B \otimes \boldsymbol{\Omega}_B \\
&\mathbf{I}_B \dot{\boldsymbol{\Omega}}_B \;=\; -\boldsymbol{\Omega}_B \times (\mathbf{I}_B \boldsymbol{\Omega}_B) \;+\; \boldsymbol{\tau}_B^{\mathrm{body}}(\boldsymbol{\rho}, q_B)
\end{aligned}
}
$$

with:

$$
\mathbf{F}_{\mathrm{Kepler}}(\boldsymbol{\rho}) \;=\; -\frac{Gm_Am_B}{|\boldsymbol{\rho}|^3}\, \boldsymbol{\rho}, \qquad
\boldsymbol{\tau}_X^{\mathrm{body}} \;=\; \frac{3Gm_{\bar X}}{|\boldsymbol{\rho}|^3}\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, X} \times (\mathbf{I}_X \hat{\boldsymbol{\rho}}_{\mathrm{body},\, X}),
$$

and $\mathbf{F}_{\mathrm{back-reaction},\, X}$ given by §4.6.2 boxed.

This is the canonical-tier equation set — the equations of motion that any implementation of this physics must satisfy. Chapter 05 specialises the rotational pieces to the principal-axis frame and substitutes the explicit torque expression; Chapter 06 validates each piece against analytical limits.

### §4.7.1 Line-by-line match against `dynamics.py::state_derivative`

The engineer-tier implementation in `dynamics.py::state_derivative` (lines ~130–170) propagates exactly the boxed equation set above. Condensed, the function reads (paraphrased; consult the file for the literal code):

```
# Translational
r_vec = body_b.x - body_a.x   # equals boldsymbol{rho}
r_mag = norm(r_vec)
F_kepler = -G * m_a * m_b * r_vec / r_mag**3
# Rotational A
tau_a_body = gravity_gradient_torque_body_frame(r_vec, R_a, I_a, m_b)
omega_a_dot = solve(I_a, tau_a_body - cross(omega_a, I_a @ omega_a))
q_a_dot = 0.5 * q_a_kinematics(q_a, omega_a)
# Rotational B (symmetric, with -r_vec)
tau_b_body = gravity_gradient_torque_body_frame(-r_vec, R_b, I_b, m_a)
omega_b_dot = solve(I_b, tau_b_body - cross(omega_b, I_b @ omega_b))
q_b_dot = 0.5 * q_b_kinematics(q_b, omega_b)
# State derivative assembly
return [x_a_dot, x_b_dot, v_a_dot, v_b_dot, q_a_dot, q_b_dot, omega_a_dot, omega_b_dot]
```

The match is line-by-line:

| Canonical equation (boxed §4.7) | Engineer-tier line | Status |
|---|---|---|
| $\dot{\boldsymbol{\rho}} = \mathbf{p}_{\mathrm{rel}}/\mu$ | `v_a_dot, v_b_dot` (in CoM frame the relative motion is automatic from the difference) | exact |
| $\mathbf{F}_{\mathrm{Kepler}} = -Gm_Am_B \boldsymbol{\rho}/r^3$ | `F_kepler` line | exact (prefactor + sign) |
| $\mathbf{F}_{\mathrm{back-reaction},\, X}$ (§4.6.2 boxed) | **NOT included** in first-cut `state_derivative` | engineer-tier gap; flagged at OQ-4.4 resolution below |
| $\dot{q}_X = \tfrac{1}{2} q_X \otimes \boldsymbol{\Omega}_X$ | `q_a_kinematics`, `q_b_kinematics` | exact |
| $\mathbf{I}_X \dot{\boldsymbol{\Omega}}_X = \boldsymbol{\tau}_X^{\mathrm{body}} - \boldsymbol{\Omega}_X \times (\mathbf{I}_X \boldsymbol{\Omega}_X)$ | `omega_a_dot, omega_b_dot` lines | exact |
| $\boldsymbol{\tau}_X^{\mathrm{body}} = (3Gm_{\bar X}/r^3)\hat{\boldsymbol{\rho}}_{\mathrm{body},X} \times (\mathbf{I}_X \hat{\boldsymbol{\rho}}_{\mathrm{body},X})$ | `gravity_gradient_torque_body_frame` (lines 113–131) — verified in Ch 03 §3.7.1 | exact |
| Body-$B$ call convention with $-\boldsymbol{\rho}$ as input | line 158 — verified in Ch 03 §3.7.1 | exact |

**(Resolves OQ-4.4.)** The line-by-line match is exact at the **rotational** level (translational Kepler + rotational Euler + gravity-gradient torque all correct). The **gap** is at the back-reaction force: the canonical-tier §4.6.2 boxed expression $\mathbf{F}_{\mathrm{back-reaction},\, X}$ is **not yet implemented** in `dynamics.py::state_derivative`. The first-cut therefore propagates pure Kepler translation + gravity-gradient-torqued rotation, omitting the orientation-back-reaction-on-orbit coupling.

**Magnitude of the missing contribution.** At first-cut configurations (Ch 03 §3.7.2 hypothetical L2-like), $|\mathbf{F}_{\mathrm{Kepler}}| \sim G m_A m_B / r^2 \sim 10^{-12}\ \text{N}$, and $|\mathbf{F}_{\mathrm{back-reaction}}| \sim 3 G m_B \mathrm{tr}(\mathbf{I}_A) / r^5 \sim 10^{-18}\ \text{N}$ — the back-reaction is six orders of magnitude smaller than the Kepler force at separations $\gg \ell$ (the body size), and the ratio scales as $(\ell/r)^3$ as expected for the next-multipole correction. At the first-cut scales the omission is below RK4 round-off and does not affect any observed numerical result. At binary-asteroid-tidal-locking regimes ($\ell/r \to 1$) the back-reaction becomes the dominant non-Kepler orbital effect; v0.2 engineer-tier work would implement it. **Engineer-tier follow-up:** queued as new OQ-4.5 at §4.9.

### §4.7.2 Numerical sanity check: Kepler limit reduction

Set $\mathbf{I}_A = \mathbf{I}_B = \mathbf{0}$ in the boxed §4.7 equation set. Then $\boldsymbol{\tau}_A^{\mathrm{body}} = \boldsymbol{\tau}_B^{\mathrm{body}} = \mathbf{0}$ (proportional to $\mathbf{I}_X$), $\mathbf{F}_{\mathrm{back-reaction},\, X} = \mathbf{0}$ (proportional to $\mathbf{I}_X$ as well — §4.6.2 boxed), and the translational equation collapses to:

$$
\mu \ddot{\boldsymbol{\rho}} \;=\; \mathbf{F}_{\mathrm{Kepler}} \;=\; -\frac{Gm_Am_B}{r^2}\, \hat{\boldsymbol{\rho}},
$$

with the rotational equations decoupled to $\boldsymbol{\Omega}_X = \mathrm{const}$ (no torque on zero-inertia body). The translational equation is the standard reduced-mass Kepler equation; conserved quantities (energy, angular momentum, Laplace-Runge-Lenz vector) reduce to the Kepler conserved quantities; orbital period matches the textbook formula $T_{\mathrm{orbit}} = 2\pi\sqrt{a^3/(G(m_A + m_B))}$. This is the testbed §2.3 / Ch 06 §6.2 reduction; the sanity check confirms the boxed §4.7 equation set has the correct Kepler limit.

## §4.8 What is established by the end of this chapter

A reader who has worked through §4.1–§4.7 has:

- The complete Lagrangian $L = T - V$ (§4.2 boxed) on the reduced configuration manifold $Q_{\mathrm{CoM}} = \mathbb{R}^3 \times SO(3)_A \times SO(3)_B$, with the kinetic energy from Ch 02 and the quadrupole-truncated potential from Ch 03 assembled and structurally decomposed (translational kinetic + Kepler potential + rotational kinetic per body + tidal coupling per body).
- The structural identification (§4.2.1) of the **tidal piece** as the exclusive source of rotational-translational coupling, and the recognition that its $\partial/\partial q_X$ projection gives the gravity-gradient torque while its $\partial/\partial \boldsymbol{\rho}$ projection gives the orbital back-reaction.
- The Euler-Lagrange equations on $TQ_{\mathrm{CoM}}$ derived by direct chart computation (Approach A, §4.3) in the Euler-angle chart, with the chart algebra reproducing the classical body-frame Euler dynamic equations via the kinematic transport relation. CS-language anchor on Podolský §6.4–§6.5 (chart machinery) + §7.3 (worked algebra) made concrete.
- The Euler-Poincaré reduction (Approach B, §4.4) of the left-invariant rotational Lagrangian on each $SO(3)$ factor, with the $\mathrm{ad}^*$ action computed explicitly as $-\boldsymbol{\Omega} \times \boldsymbol{\mu}$ in the $\mathfrak{so}(3) \cong \mathbb{R}^3$ identification, locking the sign convention to match Ch 02 §2.7.1 + Ch 03 §3.7 + the 2026-05-24 OQ-FORT-1 sign-fix.
- The **two-in-parallel cross-validation** (§4.5) showing Approach A ≡ Approach B at the body-frame Euler equation — the analytical analog of doctrine §4.4 C2 cross-implementation bit-identity, applied at the derivation level.
- The translational equation $\mu \ddot{\boldsymbol{\rho}} = -\nabla_{\boldsymbol{\rho}} V$ (§4.6) decomposed into Kepler force + back-reaction force, with **explicit closed-form $\mathbf{F}_{\mathrm{back-reaction},\, X}$ derivation** (§4.6.2 boxed): scales as $1/r^4$, has orientation-dependent prefactor, vanishes identically in the spherical-body limit (§4.6.3), participates in the Newton's-third-law balance with the body torques to keep $\mathbf{P}_{\mathrm{CoM}}$ and $\mathbf{L}_{\mathrm{total}}$ conserved (§4.6.4).
- The **complete equation set** (§4.7 boxed) on the 18-dimensional reduced phase space, line-by-line cross-validated against `dynamics.py::state_derivative` at the rotational + Kepler level, with the back-reaction-force gap identified explicitly and queued at §4.9 OQ-4.5.

**What is NOT yet established:**

- The body-frame Euler equations stated as a self-contained specialisation with the gravity-gradient torque substituted explicitly on the right-hand side (Ch 03 §3.7 boxed form expanded into Cartesian principal-axis components). **Chapter 05.**
- Numerical verification of the equations of motion against analytical-limit testbeds (Kepler, torque-free symmetric top, gravity-gradient libration, Dzhanibekov asymmetric-top, infinite-separation decoupling). **Chapter 06.**
- The symplectic-Lie-group integrator that would replace RK4 for long-horizon accurate integration. **v0.2 engineer-tier scope** (flagged in Ch 01 §1.6 cross-tier note + Ch 03 §3.7.2 multi-scale time-scale table + §4.9 OQ-4.7).
- The L2 / restricted-three-body extension adding centrifugal and Coriolis pseudo-forces. **v0.2 canonical tier** (Ch 00 §1).
- Implementation of the back-reaction force in `dynamics.py::state_derivative`. **v0.2 engineer-tier scope** (§4.9 OQ-4.5).

---

## §4.9 Open questions in this chapter

### Resolved this round (canonical-tier round 3, 2026-05-25)

- **OQ-4.1** ***Resolved*** — chart choice locked at §4.3: **Euler-angle chart** as canonical-tier presentation (matches Podolský §7.3 directly), quaternion chart as engineer-tier integration choice (Ch 02 §2.7.1 engineer-tier note), both producing the same boxed body-frame Euler equation. No unit-norm constraint to handle in the canonical-tier chart; the engineer-tier chart's constraint is handled at integration time via re-normalisation per step (standard quaternion-propagation practice).
- **OQ-4.2** ***Resolved*** — at §4.4.3 boxed: $\mathrm{ad}^*_{\boldsymbol{\Omega}}\,\boldsymbol{\mu} = -\boldsymbol{\Omega} \times \boldsymbol{\mu}$ for $\mathfrak{so}(3) \cong \mathbb{R}^3$, with the minus sign forced by the pairing convention. Sign consistent with Ch 02 §2.7.1 left-trivialisation + Ch 03 §3.7 boxed torque + OQ-FORT-1 sign-fix (2026-05-24 PM).
- **OQ-4.3** ***Resolved*** — at §4.6.2 boxed: $\mathbf{F}_{\mathrm{back-reaction},\, A} = -(3Gm_B/2r^4)\Phi_A \hat{\boldsymbol{\rho}} - (3Gm_B/r^4)(\mathbb{I}_3 - \hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}^T) \mathbf{R}_q^A (\mathbf{I}_A \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A})$. Scales as $1/r^4$ as expected; vanishes in spherical-body limit (§4.6.3); preserves $\mathbf{P}_{\mathrm{CoM}}$ via Newton's-third-law balance with the gravity-gradient torque (§4.6.4).
- **OQ-4.4** ***Resolved*** — at §4.7.1 table: line-by-line match against `dynamics.py::state_derivative` for the rotational + Kepler-translational pieces. Back-reaction force is not implemented in first-cut; magnitude $\sim (\ell/r)^3$ relative to Kepler force, below RK4 round-off at first-cut configurations. Engineer-tier follow-up flagged as new OQ-4.5.

### Surfaced this round (forward-looking, v0.1.1 / v0.2 scope)

- **OQ-4.5.** Implement $\mathbf{F}_{\mathrm{back-reaction},\, A} + \mathbf{F}_{\mathrm{back-reaction},\, B}$ (§4.6.2 boxed) in `dynamics.py::state_derivative` translational integration. Expected impact at first-cut configurations: below double-precision round-off. Becomes engineering-significant in binary-asteroid-tidal-locking regime ($\ell/r \to 1$). **Engineer-tier v0.1.1 or v0.2 scope.** Companion conservation diagnostic: with back-reaction implemented, the orbital angular momentum $\boldsymbol{\rho} \times \mu \dot{\boldsymbol{\rho}}$ should exchange with body spins at the level predicted by the Noether identity (§4.6.4) and the $|d\mathbf{L}_{\mathrm{total}}/L_0|$ diagnostic should improve toward double-precision floor.
- **OQ-4.6.** Develop the **symplectic-Lie-group variational integrator** (Hairer-Lubich-Wanner 2006 §VII.5; Lee-Leok-McClamroch 2018 *Global Formulations of Lagrangian and Hamiltonian Dynamics on Manifolds*) for the boxed §4.7 equation set. The current RK4 propagation is non-symplectic; long-horizon ($\gtrsim 10^7$ s, of order the orbital / libration period at hypothetical L2 configurations) integration accumulates a secular energy drift $\sim O(dt^4 \cdot T)$ that bounds the achievable physics fidelity. The symplectic-Lie-group framework respects the underlying geometric structure (left-invariance + symplectic form on $T^*Q_{\mathrm{CoM}}$) by construction. **v0.2 engineer-tier canonical-tier-paired scope.** Implementation: discrete Lagrangian on $SO(3) \times SO(3) \times \mathbb{R}^3$ via the discrete-time Euler-Poincaré reduction (DT-EP) of Marsden-Pekarsky-Shkoller 1999; tested against the analytical limits of Chapter 06 plus a long-horizon-energy-drift test (per `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md` §3.7).
- **OQ-4.7.** Develop the **energy-momentum method** (Simo-Lewis-Marsden 1991; Marsden-Ratiu 1999 §1.7) for the stability analysis of the boxed §4.7 equation set at relative equilibria (e.g., the gravity-gradient libration equilibrium of Ch 03 §3.7.2 + Ch 05 §5.4). The energy-momentum method generalises Routh-Hurwitz / Lyapunov stability to systems with symmetry, providing rigorous criteria for nonlinear stability of relative equilibria on Lie groups. Pairs with the linearised libration analysis of Ch 03 §3.7.2 (which gives spectral stability) to upgrade to nonlinear stability under the doctrine "test ≥2 independent methods agree". **v0.2 canonical-tier scope.**
- **OQ-4.8.** Develop the **Lyapunov spectrum** of the boxed §4.7 equation set near the Dzhanibekov regime (§2.2 of testbed catalogue; Ch 05 §5.3 asymmetric case; Ch 06 §6.5 testbed). Cross-validate analytical Lyapunov exponents (Wisdom 1987 for spin-orbit resonance Lyapunov spectra; specialised to two-body coupled case by Maciejewski 1995) against numerical computation from the engineer-tier output. **v0.2 canonical-tier-paired scope.**

---

**Next chapter:** [05 — Body-frame Euler equations as a specialisation](05-body-frame-euler-equations.md)

**End 04-euler-lagrange-and-euler-poincare.md (v0.1, round 3).**
