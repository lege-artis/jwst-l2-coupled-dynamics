# Chapter 1 — Configuration manifold and variational principle

> **Prerequisite reading.** Chapter 00 (introduction). Familiarity with elementary multivariable calculus and the notion of a smooth manifold. No prior exposure to symplectic geometry required; structures are introduced as needed.
>
> **Status.** v0.1 — first 3-5 pages of intended ~15-page chapter. Subsequent sections will be added in dedicated Opus authoring sessions.

---

## §1.1 What this chapter does

This chapter sets the geometric stage for everything that follows. We identify the **configuration space** of the two-rigid-body system — the set of all possible *positions and orientations* that the pair of bodies can occupy — and equip it with the structure needed to write Lagrangian mechanics on it: a tangent bundle, a smooth Lagrangian function, and a variational principle. By the end of the chapter, the Euler-Lagrange equations of motion for our system will be ready to derive (Chapter 04), and the question "what specific frame should we calculate in?" will be a deliberate choice rather than an assumption.

The treatment is **frame-invariant**: we describe the system in language that does not privilege any particular reference frame. When we eventually choose a frame (Chapter 05 onwards), we do so because the equations are simplest there, not because we forgot we had a choice.

## §1.2 The configuration of one rigid body

Before we describe two coupled bodies, we describe one.

A **rigid body** is a collection of mass points $\{m_\alpha\}$ at positions $\{\mathbf{r}_\alpha\}$ such that all pairwise distances $|\mathbf{r}_\alpha - \mathbf{r}_\beta|$ are constant in time. By a classical theorem (Goldstein 2002 §4.1; Arnold 1989 §3.27), the configuration of a rigid body in three-dimensional Euclidean space is completely specified by:

1. The position of any one chosen reference point (conventionally, the centre of mass) in some inertial frame: a point $\mathbf{R} \in \mathbb{R}^3$.
2. The orientation of the body relative to that inertial frame: a rotation $q \in SO(3)$.

Here $SO(3)$ denotes the **special orthogonal group** in three dimensions — the set of $3 \times 3$ real matrices $\mathbf{R}_q$ satisfying $\mathbf{R}_q^T \mathbf{R}_q = \mathbf{I}$ and $\det(\mathbf{R}_q) = +1$. Geometrically, $SO(3)$ is the set of rotations of three-space that preserve handedness. As a smooth manifold, $SO(3)$ has dimension three; one common parameterisation is by Euler angles $(\phi, \theta, \psi)$, another is by unit quaternions $q \in S^3 \subset \mathbb{R}^4$ with the identification $q \sim -q$ (double cover; see Marsden & Ratiu 1999 §9.1). The first-cut engineer-tier code uses unit quaternions; the canonical-tier derivation is parameterisation-independent.

The configuration of one rigid body is therefore an element of:

$$
SE(3) = \mathbb{R}^3 \times SO(3),
$$

the **special Euclidean group** of rigid motions in three-space. It is a six-dimensional Lie group. The first factor describes translation, the second rotation.

> **Remark.** The semidirect-product structure of $SE(3) = \mathbb{R}^3 \rtimes SO(3)$ is essential when we discuss how rotation and translation compose (rotating a translated object is not the same as translating a rotated one). We will need this structure explicitly in Chapter 04 when we discuss reduction by symmetries. For now, $SE(3) \cong \mathbb{R}^3 \times SO(3)$ as a manifold suffices.

## §1.3 The configuration of two rigid bodies

For two rigid bodies $A$ and $B$, each independently positioned and oriented in three-space, the configuration space is the Cartesian product:

$$
Q \;=\; SE(3)_A \times SE(3)_B \;=\; \big(\mathbb{R}^3 \times SO(3)\big)_A \;\times\; \big(\mathbb{R}^3 \times SO(3)\big)_B.
$$

A point in $Q$ is a 4-tuple:

$$
q \;=\; (\mathbf{R}_A,\ q_A,\ \mathbf{R}_B,\ q_B), \qquad \mathbf{R}_A, \mathbf{R}_B \in \mathbb{R}^3, \quad q_A, q_B \in SO(3).
$$

The dimension of $Q$ is $6 + 6 = 12$. This is the expected count: three positions for $A$'s centre of mass, three angles for $A$'s orientation, the same for $B$. The first-cut engineer-tier code (`dynamics.py`, `BodyState` class) carries 13 components per body (3 position + 3 velocity + 4 quaternion + 3 angular velocity = 13), which is the *phase-space* count, twice the configuration count plus the quaternion-vs-rotation-matrix overcount.

### §1.3.1 An immediate symmetry — and what we do with it

The configuration space $Q$ has a natural action of the group $SE(3)$ itself on it: if $g = (\mathbf{a}, h) \in SE(3)$ is a global rigid motion (translation by $\mathbf{a}$, rotation by $h$), then

$$
g \cdot q \;=\; \big(h \mathbf{R}_A + \mathbf{a},\ h \cdot q_A,\ h \mathbf{R}_B + \mathbf{a},\ h \cdot q_B\big)
$$

describes the same physical configuration as $q$, just viewed from a translated and rotated reference frame.

This is the formal expression of the statement that **physics does not depend on choice of inertial reference frame**. The action of $SE(3)$ on $Q$ takes a configuration to another configuration that is physically equivalent.

The quotient $Q / SE(3)$ is the space of physically distinct configurations. It is a six-dimensional manifold (twelve minus the six dimensions of the symmetry group). We will not work directly in the quotient — for computational purposes, choosing a frame to fix the symmetry is simpler — but the existence of this quotient is what lets us choose the **centre-of-mass frame** later (Chapter 02) without loss of generality. The CoM frame fixes three of the six translational/rotational symmetry parameters; the remaining three are fixed by aligning axes with the relative position vector or by chosing a body-fixed orientation. None of these choices changes the physics.

> **Convention adopted from this point.** Unless otherwise stated, we work in an **inertial frame** with origin at the centre of mass of the two-body system. Translational invariance of the equations of motion implies that this CoM is itself inertial. The total linear momentum $\mathbf{P}_{\text{total}} = m_A \dot{\mathbf{R}}_A + m_B \dot{\mathbf{R}}_B$ is constant; we set it to zero by frame choice. Six conserved quantities ($\mathbf{R}_{\text{CoM}}$ position and $\mathbf{P}_{\text{total}}$ value) are thereby trivialised. The remaining dynamics live on a six-dimensional reduced configuration space, parameterised conveniently as $(\Delta \mathbf{R},\ q_A,\ q_B)$ where $\Delta \mathbf{R} = \mathbf{R}_B - \mathbf{R}_A$ is the relative position vector and $q_A$, $q_B$ remain as before. We will use this reduction freely from Chapter 02 onward.

## §1.4 Velocities — the tangent bundle

Lagrangian mechanics requires not just positions but velocities. The set of all (position, velocity) pairs at a given configuration is the **tangent space** to the configuration manifold at that configuration; the disjoint union of tangent spaces over all configurations is the **tangent bundle** $TQ$.

For our system, the tangent bundle has dimension $2 \cdot 12 = 24$. A point in $TQ$ is a configuration $q$ together with a velocity vector:

$$
(q,\ \dot{q}) \;=\; \big( (\mathbf{R}_A, q_A, \mathbf{R}_B, q_B),\ (\dot{\mathbf{R}}_A, \dot{q}_A, \dot{\mathbf{R}}_B, \dot{q}_B) \big).
$$

The translational velocity components $\dot{\mathbf{R}}_A, \dot{\mathbf{R}}_B \in \mathbb{R}^3$ are unproblematic — they live in the tangent space of $\mathbb{R}^3$, which is $\mathbb{R}^3$ itself. The rotational velocity components $\dot{q}_A, \dot{q}_B$ require care because $SO(3)$ is not a vector space; its tangent space at a point is *not* the same as the manifold itself.

The tangent space $T_{q_A} SO(3)$ at an orientation $q_A \in SO(3)$ is a three-dimensional vector space. It can be identified with the Lie algebra $\mathfrak{so}(3)$ of antisymmetric $3 \times 3$ matrices, which in turn is naturally isomorphic to $\mathbb{R}^3$ via the **hat map**:

$$
\hat{}\colon \mathbb{R}^3 \to \mathfrak{so}(3), \qquad
\boldsymbol{\omega} = (\omega_1, \omega_2, \omega_3) \;\mapsto\;
\hat{\boldsymbol{\omega}} \;=\;
\begin{pmatrix}
0 & -\omega_3 & \omega_2 \\
\omega_3 & 0 & -\omega_1 \\
-\omega_2 & \omega_1 & 0
\end{pmatrix}.
$$

This map is an isomorphism of Lie algebras (the cross product on $\mathbb{R}^3$ becomes matrix commutator on $\mathfrak{so}(3)$): $\widehat{\mathbf{a} \times \mathbf{b}} = [\hat{\mathbf{a}}, \hat{\mathbf{b}}]$.

The physical content of the hat map: given a rotation matrix $\mathbf{R}_q(t)$ depending on time, its time derivative is $\dot{\mathbf{R}}_q = \hat{\boldsymbol{\omega}}\, \mathbf{R}_q$, where $\boldsymbol{\omega} \in \mathbb{R}^3$ is the **angular velocity in the inertial frame**. Equivalently, $\dot{\mathbf{R}}_q = \mathbf{R}_q\, \hat{\boldsymbol{\Omega}}$ where $\boldsymbol{\Omega} = \mathbf{R}_q^T \boldsymbol{\omega}$ is the **angular velocity in the body frame**. Both representations describe the same physical rotation rate; the choice between them is a question of which frame is most convenient.

> **Code cross-reference.** The first-cut `dynamics.py::q_kinematics_matrix` computes $\dot{q}$ for a unit-quaternion representation of $\mathbf{R}_q$; the formula is $\dot{q} = \tfrac{1}{2}\, q \cdot (0, \boldsymbol{\Omega})$ in quaternion algebra. This is the body-frame ($\boldsymbol{\Omega}$) version of the relation above. Derivation: take the time derivative of $\mathbf{R}_q^T \mathbf{R}_q = \mathbf{I}$ to get $\dot{\mathbf{R}}_q^T \mathbf{R}_q + \mathbf{R}_q^T \dot{\mathbf{R}}_q = 0$, hence $\mathbf{R}_q^T \dot{\mathbf{R}}_q$ is antisymmetric and equals $\hat{\boldsymbol{\Omega}}$ for some $\boldsymbol{\Omega} \in \mathbb{R}^3$. Substitute the quaternion parameterisation and reduce.

Practically, this means that **velocities on our 12-dimensional configuration manifold are described by 12 numbers**: three translational velocity components for each body's CoM, and three angular velocity components for each body's rotation (either inertial-frame or body-frame, by choice). The tangent bundle $TQ$ is twelve-dimensional in configuration plus twelve in velocity, for twenty-four total.

## §1.5 The Lagrangian and the variational principle

We have a configuration manifold $Q$ and its tangent bundle $TQ$. **Lagrangian mechanics** consists of:

1. A smooth function $L \colon TQ \times \mathbb{R} \to \mathbb{R}$ — the **Lagrangian** — assigning a real number to each (configuration, velocity, time) triple. For our autonomous (time-independent) system, $L = L(q, \dot{q})$ has no explicit time dependence.
2. A variational principle: the actual physical trajectory between two configurations $q(t_1)$ and $q(t_2)$ is the one that makes the **action functional**

$$
S[q] \;=\; \int_{t_1}^{t_2} L(q(t), \dot{q}(t))\ dt
$$

extremal under small variations of the path $q(t)$ that vanish at the endpoints. This is **Hamilton's principle** (Goldstein §2.1; Thorne & Blandford Vol I §7.2; Feynman Vol II Ch 19; Arnold §3.13). For mechanical systems with kinetic energy $T$ and potential energy $V$, the Lagrangian takes the form $L = T - V$.

> **Why a variational principle?** Three answers, in increasing order of depth:
> 1. **Operationally**: variational formulations make symmetries and conservation laws transparent via Noether's theorem (§1.6 below).
> 2. **Geometrically**: the action functional and its extremals are coordinate-independent objects; the resulting equations of motion are guaranteed to have the same form in any choice of generalised coordinates.
> 3. **Foundationally**: extending classical mechanics to quantum mechanics is more natural in the Lagrangian/action formulation (Feynman path integrals) than in the Newtonian force-and-acceleration formulation. We do not need this perspective for the canonical tier here, but its existence is a strong sanity check that this is the right starting point.

The extremality condition $\delta S = 0$ for variations $\delta q(t)$ vanishing at $t_1, t_2$ yields the **Euler-Lagrange equations**:

$$
\frac{d}{dt} \frac{\partial L}{\partial \dot{q}^i} \;-\; \frac{\partial L}{\partial q^i} \;=\; 0, \qquad i = 1, \ldots, \dim Q.
$$

Here $q^i$ are local coordinates on $Q$ and $\dot{q}^i$ their time derivatives. The derivation of this equation from $\delta S = 0$ is standard (Goldstein §2.3; Arnold §3.14); we omit it because no problem-specific structure of our two-body system enters.

For our specific system, Chapters 02 (kinetic energy on $SE(3) \times SE(3)$), 03 (mutual gravitational potential), and 04 (assembling $L = T - V$ and applying the Euler-Lagrange operator) carry out the construction in detail.

## §1.6 Noether's theorem — symmetries become conservation laws

A continuous symmetry of the Lagrangian — that is, a one-parameter family of transformations $q \mapsto q_\epsilon$ such that $L(q_\epsilon, \dot{q}_\epsilon) = L(q, \dot{q})$ for all $\epsilon$ near zero — implies a conserved quantity along physical trajectories. This is **Noether's theorem** (Goldstein §13.7; Arnold §A.7; Marsden & Ratiu §11.2).

For our two-body system, the immediate symmetries of the Lagrangian (which we will verify formally in Chapter 04) are:

| Symmetry | Conserved quantity | Conservation diagnostic in engineer-tier |
|---|---|---|
| Time translation $t \to t + \epsilon$ | Total energy $E = T + V$ | `dynamics.py::total_kinetic_energy + total_potential_energy`; first-cut achieves $\|dE/E_0\| = 6.2 \times 10^{-12}$ |
| Spatial translation $\mathbf{R}_A, \mathbf{R}_B \to \mathbf{R}_A + \boldsymbol{a}, \mathbf{R}_B + \boldsymbol{a}$ | Total linear momentum $\mathbf{P}$ | First-cut: $\mathbf{P}$ remains zero to machine precision since CoM frame is enforced at initialisation |
| Spatial rotation by $h \in SO(3)$ acting on both bodies | Total angular momentum $\mathbf{L}$ | `dynamics.py::total_angular_momentum`; first-cut achieves $\|d\mathbf{L}/\mathbf{L}_0\| = 2.1 \times 10^{-10}$ |

The third row is structurally important: total angular momentum is conserved precisely because the Lagrangian is invariant under simultaneous global rotation of both bodies. This is *not* the same as the angular momentum of *each individual body* being conserved (each body's angular momentum exchanges with the other's via the gravity-gradient torque). The symmetry that gives total $\mathbf{L}$ conservation is the rotational invariance of the *mutual* potential energy.

> **Cross-tier note.** Noether's theorem is what *promises* the first-cut conservation diagnostics will work to machine precision in the absence of numerical error. The variational structure of the underlying equations of motion guarantees this. RK4 (the engineer-tier integrator) is non-symplectic and therefore introduces a small per-step error in energy and angular-momentum conservation; the observed $10^{-12}$ / $10^{-10}$ scales are the cumulative integration error over 600 seconds. A symplectic Lie-group integrator (Hairer-Lubich-Wanner 2006 §VII.5) would push energy conservation closer to round-off; we discuss this in Chapter 06.

## §1.7 What is established by the end of this chapter

A reader who has worked through §1.1–§1.6 has:

- The configuration space $Q = SE(3) \times SE(3)$ of our system, with explicit dimension 12.
- The tangent bundle $TQ$ of dimension 24, with velocities parameterised by translational velocities $\dot{\mathbf{R}}_A, \dot{\mathbf{R}}_B$ and angular velocities $\boldsymbol{\omega}_A, \boldsymbol{\omega}_B$ (either in inertial or body frame, our choice).
- The variational principle that selects the physical trajectory between given endpoints: $\delta S = 0$ for $S = \int L\, dt$.
- The Euler-Lagrange equations as the local form of the variational principle.
- Noether's theorem connecting symmetries of $L$ to conservation laws, ready to apply once $L$ is specified.

**What is NOT yet established:**

- The kinetic energy $T$ explicitly as a function on $TQ$. This requires the inertia tensors of the two bodies. **Chapter 02.**
- The potential energy $V$ explicitly. This requires the mutual gravitational potential as a function of the configuration. **Chapter 03.**
- The Euler-Lagrange equations in explicit coordinate form for our specific Lagrangian. **Chapter 04.**
- The Euler equations in body frame as a specialisation by Euler-Poincaré reduction. **Chapter 05.**
- Validation against analytical limits. **Chapter 06.**

> **Reference anchor.** Configuration manifolds, tangent bundles, and the variational principle: Arnold 1989 *Mathematical Methods of Classical Mechanics* §§1.4, 4.A (configuration manifolds + Lagrangian systems on $TM$); Marsden & Ratiu 1999 *Introduction to Mechanics and Symmetry* §§1, 7, 9.1 ($SO(3)$ and $SE(3)$ as Lie groups, with tangent-bundle structure); Sussman & Wisdom 2015 *Structure and Interpretation of Classical Mechanics* §1 (computational rendering of the same material). For $SO(3)$ as a Lie group specifically, Holm-Schmah-Stoica 2009 §6 (preparation for Chapter 04's Euler-Poincaré reduction).
>
> **Supplementary anchor for Czech-reading audiences.** Podolský 2018 Kapitola 6 *Kinematika tuhého tělesa* provides the elementary coordinate machinery in CS: §6.1 *Vektory a tenzory* (tensor algebra in Euclidean space, classical vs modern definitions), §6.2 *Relativita otáčivého pohybu* (inertial vs body-fixed basis $\mathbf{e}'_i = A_{ik}(t)\mathbf{e}_k$ with orthogonal $A$), §6.4 *Eulerovy úhly* (intrinsic-ZXZ parametrisation $A = BCD$ with Eqs. (6.14)–(6.16) for the three elementary rotation matrices). Podolský does *not* abstract to the Lie-group level — he treats rotations as orthogonal $3 \times 3$ matrices parameterised by Euler angles or unit quaternions, not as elements of the smooth manifold $SO(3)$ with its tangent-bundle structure and left-invariant vector fields. The Lie-group abstraction used in this chapter (and exploited in Chapter 02 §2.7.1 and Chapter 04 §4.4) is *modern Marsden-Ratiu material*, not in Podolský. For a CS-reading audience, Podolský §6.2–§6.4 is the prerequisite-coordinate reading; the Lie-group abstraction here is the modern reformulation that connects the elementary coordinate machinery to the Euler-Poincaré reduction that Chapter 04 will apply.

---

**Section 1.8 — Open questions in this chapter** (placeholders, to be resolved in subsequent authoring passes):

- **OQ-1.1.** Make the global $SE(3)$ action on $Q$ from §1.3.1 fully explicit, including the formula for the action on quaternion-parameterised orientations $q_A, q_B$.
- **OQ-1.2.** Give the formal definition of $TQ$ as a smooth manifold, including charts in a neighbourhood of an arbitrary point. (Currently described informally.)
- **OQ-1.3.** Provide the derivation of $\dot{\mathbf{R}}_q = \hat{\boldsymbol{\omega}}\, \mathbf{R}_q$ from first principles, including the connection between left-invariant and right-invariant vector fields on $SO(3)$. Reference: Marsden & Ratiu §9.

---

**Next chapter:** [02 — Kinetic energy on $SE(3) \times SE(3)$ in intrinsic form](02-kinetic-energy.md) *(queued)*

**End 01-configuration-manifold.md**
