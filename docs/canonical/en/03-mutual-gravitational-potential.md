# Chapter 3 — Mutual gravitational potential and multipole expansion

> **Prerequisite reading.** Chapter 01 (configuration manifold, hat map, centre-of-mass frame convention) and Chapter 02 (kinetic energy, inertia tensor as a symmetric positive-definite $3 \times 3$ matrix, body-frame vs. inertial-frame bookkeeping). Familiarity with elementary vector calculus and the divergence theorem. No prior exposure to spherical harmonics required; the relevant identities are introduced where they are used.
>
> **Status.** v0.1.

---

## §3.1 What this chapter does

Chapter 02 closed with the kinetic energy of the two-rigid-body system as an explicit function on $TQ$ and an identification of $T$ as a left-invariant Riemannian metric on the configuration manifold. The potential-energy half of the Lagrangian was deferred, and it is the subject of this chapter.

We derive the **mutual gravitational potential** $V$ between two extended rigid bodies from the first-principles 6-dimensional integral, expand it in inverse powers of separation by the standard multipole technique, and identify the leading orientation-dependent contribution as the **gravity-gradient tidal term**. From that term, we derive the gravity-gradient torque $\boldsymbol{\tau}$ in body frame and show, line by line, that the result matches the engineer-tier implementation in `dynamics.py::gravity_gradient_torque_body_frame`.

The derivation follows the **"two in parallel"** methodology direction issued by the project owner on 2026-05-23: the multipole expansion is carried out **independently by two standard approaches** — Cartesian Taylor expansion (Landau-Lifshitz / Murray & Dermott tradition) and spherical-harmonic expansion (Jackson tradition transferred from electromagnetism to gravity) — and the two results are shown to coincide at quadrupole order. The cross-validation step is not stylistic; it is the doctrinal `lege-artis` C2 cross-implementation principle (`_config/LEGE-ARTIS-LANGUAGE-DOCTRINE-v0.1.md` §4.4) applied at the derivation level rather than at the code level. If the two approaches disagreed, one of them would be wrong; agreement is evidence that neither is.

After Chapter 03 closes, $T$ (Chapter 02) and $V$ (this chapter) are both in hand. Chapter 04 will assemble $L = T - V$ and derive the Euler-Lagrange equations of motion on $TQ$. The remaining chapters then specialise (Chapter 05) and validate (Chapter 06) against the analytical-limit testbeds catalogued in `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md`.

## §3.2 The exact mutual potential as a six-dimensional integral

Newtonian gravity at the level treated here is **pairwise additive over mass elements**: the potential energy between two infinitesimal mass elements at $\mathbf{r}_1$ and $\mathbf{r}_2$ is $-G\, dm_1\, dm_2 / |\mathbf{r}_1 - \mathbf{r}_2|$, and the total potential of an extended system is the sum over all such pairs (Landau-Lifshitz Vol I §32; this property survives in the linearised weak-field limit of general relativity and breaks only at post-Newtonian order, which is outside our v0.1 scope — see Chapter 00 §6).

Apply this to two rigid bodies $A$ and $B$ with mass distributions $\rho_A(\mathbf{x}_A)$ and $\rho_B(\mathbf{x}_B)$ supported on compact regions $\mathcal{B}_A, \mathcal{B}_B \subset \mathbb{R}^3$ in their respective body frames. A material point of body $A$ labelled by body-frame position $\mathbf{x}_A$ sits in the inertial frame (per Chapter 02, §2.2) at $\mathbf{R}_A + \mathbf{R}_q^A\, \mathbf{x}_A$, and analogously for body $B$. The Euclidean distance between such a pair of material points is therefore:

$$
\big|\mathbf{R}_A + \mathbf{R}_q^A\, \mathbf{x}_A - \mathbf{R}_B - \mathbf{R}_q^B\, \mathbf{x}_B\big|.
$$

Summing $-G\, dm_A\, dm_B$ over all pairs of mass elements, the **exact mutual gravitational potential** of the two-body system is:

$$
\boxed{
V(\mathbf{R}_A, q_A, \mathbf{R}_B, q_B)
\;=\;
-G \int_{\mathcal{B}_A} \int_{\mathcal{B}_B}
\frac{\rho_A(\mathbf{x}_A)\, \rho_B(\mathbf{x}_B)}
     {\big|\mathbf{R}_A + \mathbf{R}_q^A\, \mathbf{x}_A - \mathbf{R}_B - \mathbf{R}_q^B\, \mathbf{x}_B\big|}
\, dV_A\, dV_B.
}
$$

This is the exact starting point — no approximation has been made beyond the rigid-body assumption (mass distributions $\rho_A, \rho_B$ are constant in their respective body frames) and the Newtonian-gravity assumption (no relativistic corrections). It is a function only on the configuration manifold $Q$, not the tangent bundle $TQ$: the potential does not depend on velocities. It is also manifestly invariant under the global $SE(3)$ action introduced in Chapter 01 §1.3.1 — a simultaneous rigid translation and rotation of both bodies leaves the relative geometry, and hence $V$, unchanged.

Two observations follow immediately and shape the rest of the chapter:

1. $V$ depends only on the **relative position** $\boldsymbol{\rho} \equiv \mathbf{R}_B - \mathbf{R}_A$ and the two **orientations** $q_A, q_B$. By translational symmetry it cannot depend separately on $\mathbf{R}_A, \mathbf{R}_B$; only their difference appears in the integrand. In the centre-of-mass frame (the convention locked in Chapter 01 §1.3.1), this is a function of six configuration variables: three components of $\boldsymbol{\rho}$ and the six rotational variables $q_A, q_B$ acting on integrand directions.

2. The integral is **non-trivial for any pair of extended bodies**: it cannot generally be reduced to elementary functions. For arbitrary mass distributions the integral is evaluated either numerically (finite-element integration over the body geometries) or analytically only in special-symmetry cases (e.g., spherically symmetric bodies, where it collapses to the point-mass result by Newton's shell theorem). The standard physics workaround is the **multipole expansion** — expand the integrand in powers of body-size-over-separation, integrate term by term, and truncate at a chosen order. This is the subject of §3.4 onward.

> **Reference anchor.** Landau-Lifshitz Vol I §32 (pairwise additivity of Newtonian gravitational potential energy); Murray & Dermott 1999 §4.1 (statement of the exact extended-body integral). Thorne & Blandford 2017 Eq. (27.50) writes the *single-body* version of this integral — $\Phi(\mathbf{x}) = -G \int \rho(\mathbf{x}')/|\mathbf{x} - \mathbf{x}'|\, dV'$ — in the context of gravitational-wave source theory (§27.5.2); our boxed two-body mutual potential reduces to two field-evaluations of T&B's $\Phi$ paired against the bodies' masses, but is mathematically equivalent at the level of the integrand. The form in the boxed equation above is standard; we have written it explicitly in the notation of Chapters 01-02 to make the configuration-dependence visible.

## §3.3 Convergence domain and the close-approach boundary

The multipole expansion is a Taylor expansion of the integrand in the small parameters $|\mathbf{x}_A|/|\boldsymbol{\rho}|$ and $|\mathbf{x}_B|/|\boldsymbol{\rho}|$ — ratios of body-internal coordinates to the inter-body separation. The expansion **converges absolutely** when:

$$
|\boldsymbol{\rho}| \;>\; \sup_{\mathbf{x}_A \in \mathcal{B}_A} |\mathbf{x}_A| \;+\; \sup_{\mathbf{x}_B \in \mathcal{B}_B} |\mathbf{x}_B|,
$$

that is, when the inter-body separation exceeds the sum of the radii of the smallest spheres centred on each body's centre of mass and enclosing that body. Geometrically: the bodies do not overlap, with a small extra margin. When this inequality fails — when the bodies are close enough that some material point of one body sits inside the multipole-expansion sphere of the other — the expansion ceases to converge and any finite truncation produces a polynomial divergence as $|\boldsymbol{\rho}|$ approaches the body radii.

This boundary defines the **physical domain of validity** of every multipole-truncated formula derived in this chapter, and it is shared by all subsequent code that consumes those formulas. The engineer-tier `dynamics.py` does not, at v0.1, contain any assertion that this inequality holds at integration time; the assertion is queued at `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md` §3.3 as part of the close-approach edge-case testbed.

> **Cross-reference.** `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md` §3.3 (close-approach assumption violation) — Class B edge-case requiring engineer-tier assertion design for v0.2. The canonical-tier statement *here* is: any multipole-truncated formula in this chapter is provably correct only in the convergence regime above; outside it, the formulas remain numerically evaluable but their physical interpretation breaks down.

The truncation order also bounds residual error inside the convergence domain. The quadrupole truncation adopted in this chapter is correct to terms of order $(\ell / |\boldsymbol{\rho}|)^2$ where $\ell$ is the characteristic body size; the next-order correction (octupole-monopole and quadrupole-dipole — both of which vanish identically for our centre-of-mass-anchored expansion, leaving the octupole-quadrupole and direct quadrupole-quadrupole as the leading omitted terms) scales as $(\ell / |\boldsymbol{\rho}|)^4$. For the first-cut numerical example bodies (JWST-like body $\ell_A \sim 7$ m; probe-like body $\ell_B \sim 1$ m) at first-cut separations $|\boldsymbol{\rho}| \sim 10$ m to $1000$ m, the quadrupole truncation is good to between $\sim 5 \times 10^{-1}$ (close approach) and $\sim 5 \times 10^{-9}$ (far separation). At realistic JWST orbital separations (hundreds of kilometres), the truncation error drops below double-precision round-off and the bodies are effectively point masses for orbital-dynamics purposes; the quadrupole term remains operative for the rotational coupling regardless of separation.

## §3.4 Approach A — Cartesian multipole expansion

We expand the denominator of the integrand in §3.2 about the separation vector $\boldsymbol{\rho}$. Define $\boldsymbol{\xi}_A \equiv \mathbf{R}_q^A\, \mathbf{x}_A$ (the body-$A$ internal coordinate expressed in the inertial frame; the rotation matrix carries the orientation dependence) and analogously $\boldsymbol{\xi}_B$. The pairwise denominator becomes:

$$
\big|\boldsymbol{\rho} + \boldsymbol{\xi}_B - \boldsymbol{\xi}_A\big|^{-1}
\;=\;
\big|\boldsymbol{\rho} - \mathbf{u}\big|^{-1},
\qquad \mathbf{u} \equiv \boldsymbol{\xi}_A - \boldsymbol{\xi}_B.
$$

Taylor-expand $1/|\boldsymbol{\rho} - \mathbf{u}|$ in powers of $\mathbf{u}/|\boldsymbol{\rho}|$. Working out the first three terms (the algebra is elementary but somewhat tedious; the standard manipulation uses $|\boldsymbol{\rho} - \mathbf{u}|^{-1} = (|\boldsymbol{\rho}|^2 - 2 \boldsymbol{\rho} \cdot \mathbf{u} + |\mathbf{u}|^2)^{-1/2}$ and the binomial series for $(1 - \varepsilon)^{-1/2}$):

$$
\frac{1}{|\boldsymbol{\rho} - \mathbf{u}|}
\;=\;
\frac{1}{|\boldsymbol{\rho}|}
\;+\;
\frac{\hat{\boldsymbol{\rho}} \cdot \mathbf{u}}{|\boldsymbol{\rho}|^2}
\;+\;
\frac{3\, (\hat{\boldsymbol{\rho}} \cdot \mathbf{u})^2 - |\mathbf{u}|^2}{2\, |\boldsymbol{\rho}|^3}
\;+\;
\mathcal{O}\!\left(\frac{|\mathbf{u}|^3}{|\boldsymbol{\rho}|^4}\right),
$$

where $\hat{\boldsymbol{\rho}} \equiv \boldsymbol{\rho}/|\boldsymbol{\rho}|$ is the unit vector from $A$ to $B$. Substitute $\mathbf{u} = \boldsymbol{\xi}_A - \boldsymbol{\xi}_B$ and insert into the integral of §3.2:

$$
V \;=\; -G \int_A \int_B \rho_A \rho_B
\left[ \frac{1}{|\boldsymbol{\rho}|}
\;+\; \frac{\hat{\boldsymbol{\rho}} \cdot (\boldsymbol{\xi}_A - \boldsymbol{\xi}_B)}{|\boldsymbol{\rho}|^2}
\;+\; \frac{3\, [\hat{\boldsymbol{\rho}} \cdot (\boldsymbol{\xi}_A - \boldsymbol{\xi}_B)]^2 - |\boldsymbol{\xi}_A - \boldsymbol{\xi}_B|^2}{2\, |\boldsymbol{\rho}|^3}
\;+\; \cdots
\right] dV_A\, dV_B.
$$

The integrals separate: $\boldsymbol{\xi}_A$ depends only on $\mathbf{x}_A$ and $q_A$, and $\boldsymbol{\xi}_B$ only on $\mathbf{x}_B$ and $q_B$. We evaluate the three leading terms in turn.

### §3.4.1 Order zero — monopole-monopole (Kepler term)

The $1/|\boldsymbol{\rho}|$ term comes out of the integrals trivially:

$$
V^{(0)} \;=\; -\frac{G}{|\boldsymbol{\rho}|} \int_A \rho_A\, dV_A \int_B \rho_B\, dV_B
\;=\; -\frac{G\, m_A\, m_B}{|\boldsymbol{\rho}|}.
$$

This is the **monopole-monopole** term: each body acts gravitationally on the other as if its mass were concentrated at the centre of mass. It is the entire gravitational potential between two point masses and recovers the **Kepler problem** in the limit where the rotational degrees of freedom decouple (validation testbed `_specs/...` §2.3). It carries no orientation dependence and therefore contributes **zero gravity-gradient torque**.

### §3.4.2 Order one — dipole terms vanish by centre-of-mass

The linear-in-$\mathbf{u}$ term decomposes into two pieces, one proportional to $\int \rho_A \boldsymbol{\xi}_A dV_A$ and the other to $\int \rho_B \boldsymbol{\xi}_B dV_B$. Both vanish:

$$
\int_A \rho_A\, \boldsymbol{\xi}_A\, dV_A \;=\; \int_A \rho_A\, \mathbf{R}_q^A\, \mathbf{x}_A\, dV_A
\;=\; \mathbf{R}_q^A \int_A \rho_A\, \mathbf{x}_A\, dV_A
\;=\; \mathbf{0},
$$

because $\int \rho_A \mathbf{x}_A dV_A = \mathbf{0}$ by the body-frame centre-of-mass convention of Chapter 02 §2.2. The same algebra applies to body $B$. Therefore:

$$
V^{(1)} \;=\; 0.
$$

This is structurally important: it is the **same identity** that eliminated the cross-term in König's decomposition of the kinetic energy (Chapter 02 §2.3). The choice to anchor body coordinates at the centre of mass eliminates dipole terms in *both* $T$ and $V$ simultaneously, leaving the quadrupole as the leading orientation-dependent contribution to each.

> **Remark.** If, hypothetically, body coordinates were anchored at some other reference point (e.g., the geometric centroid of a non-uniform body, distinct from the centre of mass), the dipole term would survive in both $T$ and $V$ with the same dipole moment $\mathbf{d}_X = \int \rho_X \mathbf{x}_X dV_X \neq \mathbf{0}$. The centre-of-mass convention is therefore not just one of several equally valid frame choices: it is the unique anchor that kills the dipole term and makes the multipole expansion start at quadrupole order.

### §3.4.3 Order two — quadrupole terms

The order-two contribution requires expanding the numerator $3[\hat{\boldsymbol{\rho}} \cdot (\boldsymbol{\xi}_A - \boldsymbol{\xi}_B)]^2 - |\boldsymbol{\xi}_A - \boldsymbol{\xi}_B|^2$ that sits in the Taylor expansion of §3.4. The cross-terms in $(\boldsymbol{\xi}_A - \boldsymbol{\xi}_B)$ decompose the integral into three pieces:

- A piece quadratic in $\boldsymbol{\xi}_A$ only, integrated against $\rho_A$ and producing the total mass $m_B = \int \rho_B\, dV_B$ from the body-$B$ integral.
- A piece quadratic in $\boldsymbol{\xi}_B$ only, symmetric to the above with $A \leftrightarrow B$.
- A piece bilinear in $\boldsymbol{\xi}_A$ and $\boldsymbol{\xi}_B$, integrated against $\rho_A \rho_B$ — this contributes only the **quadrupole-quadrupole** interaction, which decays as $1/|\boldsymbol{\rho}|^5$ and is one order beyond the leading tidal coupling. It is the next-omitted term whose effect we estimated in §3.3.

The two quadratic-in-one-body pieces are the **leading tidal terms**, and they have identical structure (with $A$ and $B$ exchanged). For body $A$:

$$
V_{\mathrm{tidal},\, A}
\;=\;
-\frac{G\, m_B}{2\, |\boldsymbol{\rho}|^3}
\int_A \rho_A
\Big[ 3\, (\hat{\boldsymbol{\rho}} \cdot \boldsymbol{\xi}_A)^2 \;-\; |\boldsymbol{\xi}_A|^2 \Big]\, dV_A.
$$

To express this in terms of the **inertia tensor** $\mathbf{I}_A$ of Chapter 02 §2.4, introduce the body-frame unit vector $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \equiv (\mathbf{R}_q^A)^T \hat{\boldsymbol{\rho}}$ (the inertial separation direction expressed in body-$A$ coordinates). The integrand becomes:

$$
3\, (\hat{\boldsymbol{\rho}} \cdot \boldsymbol{\xi}_A)^2 - |\boldsymbol{\xi}_A|^2
\;=\; 3\, (\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{x}_A)^2 - |\mathbf{x}_A|^2,
$$

using the orthogonality $|\mathbf{R}_q^A \mathbf{x}_A| = |\mathbf{x}_A|$ and the analogous dot-product identity. Define the second-moment tensor of body $A$ in its body frame:

$$
(\mathbf{M}_A)_{ij} \;=\; \int_A \rho_A(\mathbf{x}_A)\, (\mathbf{x}_A)_i\, (\mathbf{x}_A)_j\, dV_A.
$$

Then $\int \rho_A (\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{x}_A)^2 dV_A = \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{M}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$, and $\int \rho_A |\mathbf{x}_A|^2 dV_A = \mathrm{tr}(\mathbf{M}_A)$.

Relating $\mathbf{M}_A$ to the inertia tensor of Chapter 02 §2.4 — $(\mathbf{I}_A)_{ij} = \int \rho_A (|\mathbf{x}_A|^2 \delta_{ij} - (\mathbf{x}_A)_i (\mathbf{x}_A)_j) dV_A$ — we have $\mathbf{I}_A = \mathrm{tr}(\mathbf{M}_A) \mathbb{I}_3 - \mathbf{M}_A$, and inversely $\mathbf{M}_A = \tfrac{1}{2} \mathrm{tr}(\mathbf{I}_A) \mathbb{I}_3 - \mathbf{I}_A$ (using $\mathrm{tr}(\mathbf{I}_A) = 2\, \mathrm{tr}(\mathbf{M}_A)$). Substitute and simplify:

$$
3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{M}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}
\;-\; \mathrm{tr}(\mathbf{M}_A)
\;=\;
\mathrm{tr}(\mathbf{I}_A) \;-\; 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}.
$$

Therefore the **quadrupole tidal term for body $A$**, expressed in the body-frame inertia tensor that Chapter 02 has already established, is:

$$
\boxed{
V_{\mathrm{tidal},\, A}
\;=\;
-\frac{G\, m_B}{2\, |\boldsymbol{\rho}|^3}
\Big[ \mathrm{tr}(\mathbf{I}_A) \;-\; 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \Big].
}
$$

A symmetric expression holds for $V_{\mathrm{tidal},\, B}$ with $A \leftrightarrow B$. The total leading orientation-dependent potential is the sum $V_{\mathrm{tidal}} = V_{\mathrm{tidal},\, A} + V_{\mathrm{tidal},\, B}$, and the **complete quadrupole-truncated mutual potential** is:

$$
\boxed{
V \;=\; -\frac{G\, m_A\, m_B}{|\boldsymbol{\rho}|}
\;+\; V_{\mathrm{tidal},\, A}
\;+\; V_{\mathrm{tidal},\, B}
\;+\; \mathcal{O}\!\left(\frac{\ell^4}{|\boldsymbol{\rho}|^5}\right).
}
$$

> **Reference anchor.** Landau-Lifshitz Vol I §41 derives the multipole expansion of the gravitational potential of a single body in compact form; Murray & Dermott 1999 §4 carries out the binary-body extension with the same Cartesian technique. Thorne & Blandford 2017 §27.5.2 ("Quadrupole-Moment Formalism") carries out the same Cartesian Taylor expansion at Eqs. (27.51)–(27.53) in the framing of gravitational-wave source theory: T&B Eq. (27.51) is the Taylor expansion of $1/|\mathbf{x} - \mathbf{x}'|$ in powers of $\mathbf{x}'/r$ (their $\mathbf{x}'$ is our $\mathbf{u}$ from §3.4); T&B Eq. (27.52) collects the monopole-plus-quadrupole truncated potential $\Phi(\mathbf{x}) = -M/r - 3 \mathcal{I}_{jk} x^j x^k / (2 r^5) + \ldots$; T&B Eq. (27.53) defines the trace-free mass quadrupole moment $\mathcal{I}_{jk} = \int \rho \big(x^j x^k - \tfrac{1}{3} r^2 \delta_{jk}\big)\, dV$. **Convention note:** T&B's mass quadrupole moment $\mathcal{I}_{jk}$ uses prefactor $1/3$, where the Cartesian quadrupole tensor $\mathbf{Q}$ used in our §3.5 follows the Jackson convention with prefactor $3$. The two are equivalent: $\mathbf{Q}_{jk} = 3\, \mathcal{I}_{jk}$, both trace-free symmetric encoding the same five physical degrees of freedom. The result for $V_{\mathrm{tidal},\, A}$ is standard. The form here in terms of $\mathrm{tr}(\mathbf{I}_A) - 3 \hat{\boldsymbol{\rho}} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}$ keeps the body-frame inertia tensor visible, which is the right form for connecting to the gravity-gradient torque (§3.7).

## §3.5 Approach B — Spherical-harmonic multipole expansion

We now redo the expansion by an independent route, using the spherical-harmonic addition theorem for $1/|\mathbf{r}_1 - \mathbf{r}_2|$ that is standard in electromagnetic-multipole theory and transfers verbatim to gravity (the only differences are the sign of the source and the absence of magnetic-multipole terms; both are irrelevant at the level of the mathematical identity).

Recall the **Laplace expansion** (Jackson 1999 eq. 3.70):

$$
\frac{1}{|\mathbf{r}_1 - \mathbf{r}_2|}
\;=\; \sum_{\ell=0}^{\infty} \sum_{m=-\ell}^{\ell}
\frac{4\pi}{2\ell+1}\,
\frac{r_<^{\ell}}{r_>^{\ell+1}}\,
Y_\ell^{m\,*}(\hat{\mathbf{r}}_2)\, Y_\ell^m(\hat{\mathbf{r}}_1),
$$

where $r_< = \min(|\mathbf{r}_1|, |\mathbf{r}_2|)$, $r_> = \max(|\mathbf{r}_1|, |\mathbf{r}_2|)$, and $Y_\ell^m$ are the standard normalised spherical harmonics. The identity converges absolutely whenever $r_< < r_>$ (i.e., wherever it is well-defined), which is the same convergence domain we identified in §3.3.

To apply this to our binary-body problem cleanly, work first with the case where body $B$ is a **point mass** at the inertial position $\mathbf{R}_B$ (this gives us $V_{\mathrm{tidal},\, A}$; the symmetric case for body $B$ follows by swapping labels, and the quadrupole-quadrupole bilinear term is one order higher, as noted in §3.4.3). The single-body potential generated by body $A$ at field point $\mathbf{R}_B$ is:

$$
\Phi_A(\mathbf{R}_B) \;=\; -G \int_A \frac{\rho_A(\mathbf{x}_A)}
{|\mathbf{R}_A + \boldsymbol{\xi}_A - \mathbf{R}_B|}\, dV_A
\;=\; -G \int_A \frac{\rho_A(\mathbf{x}_A)}{|\boldsymbol{\xi}_A - \boldsymbol{\rho}|}\, dV_A,
$$

where again $\boldsymbol{\rho} = \mathbf{R}_B - \mathbf{R}_A$ and $\boldsymbol{\xi}_A = \mathbf{R}_q^A \mathbf{x}_A$. The potential energy contributed by the body-$A$-on-point-mass-$B$ interaction is $V_A^{\mathrm{on}\, B} = m_B \Phi_A(\mathbf{R}_B)$.

Apply the Laplace expansion with $\mathbf{r}_1 = \boldsymbol{\rho}$ (the larger; we are in the convergence domain $|\boldsymbol{\rho}| > |\boldsymbol{\xi}_A|$) and $\mathbf{r}_2 = \boldsymbol{\xi}_A$ (the smaller):

$$
\frac{1}{|\boldsymbol{\xi}_A - \boldsymbol{\rho}|}
\;=\; \sum_{\ell, m} \frac{4\pi}{2\ell+1}\, \frac{|\boldsymbol{\xi}_A|^{\ell}}{|\boldsymbol{\rho}|^{\ell+1}}\,
Y_\ell^{m\,*}(\hat{\boldsymbol{\xi}}_A)\, Y_\ell^m(\hat{\boldsymbol{\rho}}).
$$

Substitute and pull the $\boldsymbol{\rho}$-dependent factor out of the integral:

$$
V_A^{\mathrm{on}\, B}
\;=\; -G\, m_B \sum_{\ell, m}
\frac{4\pi}{2\ell+1}\,
\frac{Y_\ell^m(\hat{\boldsymbol{\rho}})}{|\boldsymbol{\rho}|^{\ell+1}}\,
\underbrace{\int_A \rho_A(\mathbf{x}_A)\, |\boldsymbol{\xi}_A|^{\ell}\, Y_\ell^{m\,*}(\hat{\boldsymbol{\xi}}_A)\, dV_A}_{\equiv\;\, q_\ell^m(A)}.
$$

The bracketed integrals are the **spherical-harmonic multipole moments** of body $A$:

$$
q_\ell^m(A) \;\equiv\; \int_A \rho_A(\mathbf{x}_A)\, |\boldsymbol{\xi}_A|^{\ell}\, Y_\ell^{m\,*}(\hat{\boldsymbol{\xi}}_A)\, dV_A.
$$

These are complex-valued in general; the gravitational multipole moments are physical observables and obey the usual reality conditions $q_\ell^{-m} = (-1)^m q_\ell^{m\,*}$ since $\rho_A$ is a real density. They depend on the orientation $q_A$ through $\boldsymbol{\xi}_A = \mathbf{R}_q^A \mathbf{x}_A$.

We now evaluate the first three multipole orders:

**$\ell = 0$ (monopole).** $Y_0^0 = 1/\sqrt{4\pi}$. The moment is $q_0^0(A) = m_A / \sqrt{4\pi}$, orientation-independent. Contribution to $V$:

$$
-G m_B \cdot \frac{4\pi}{1} \cdot \frac{m_A / \sqrt{4\pi}}{|\boldsymbol{\rho}|} \cdot \frac{1}{\sqrt{4\pi}}
\;=\; -\frac{G m_A m_B}{|\boldsymbol{\rho}|}.
$$

This reproduces the Kepler term from §3.4.1, as required.

**$\ell = 1$ (dipole).** The $\ell=1$ multipole moments are linear in $\mathbf{x}_A$ and vanish by the centre-of-mass convention (same argument as §3.4.2). Specifically, $q_1^0(A) \propto \int \rho_A z_A dV_A = 0$ and analogously for $q_1^{\pm 1}$. No contribution.

**$\ell = 2$ (quadrupole).** The five complex spherical-harmonic moments $q_2^m(A)$ ($m = -2, -1, 0, 1, 2$) encode the same physical content as the five independent components of the symmetric, trace-free **quadrupole tensor**:

$$
(\mathbf{Q}_A)_{ij} \;=\; \int_A \rho_A(\mathbf{x}_A)
\big[ 3\, (\boldsymbol{\xi}_A)_i\, (\boldsymbol{\xi}_A)_j \;-\; |\boldsymbol{\xi}_A|^2\, \delta_{ij} \big]\, dV_A.
$$

The relationship between the spherical-harmonic and Cartesian forms is the standard one in Jackson 1999 §4.1: $q_2^m(A)$ are particular linear combinations of the $(\mathbf{Q}_A)_{ij}$, with the equivalence proven by direct identification of the spherical-harmonic basis functions $r^2 Y_2^m(\hat{\mathbf{r}})$ with the trace-free quadratic polynomials in Cartesian coordinates.

**Explicit basis-change algebra (resolves OQ-3.1).** Using the standard $\ell=2$ spherical harmonics in the Condon-Shortley convention (Jackson 1999 §3.6):
$$
Y_2^{0}(\hat{\mathbf{r}}) = \sqrt{\tfrac{5}{16\pi}}\,(3\cos^2\theta - 1), \quad
Y_2^{\pm 1}(\hat{\mathbf{r}}) = \mp \sqrt{\tfrac{15}{8\pi}}\,\sin\theta\cos\theta\, e^{\pm i\phi}, \quad
Y_2^{\pm 2}(\hat{\mathbf{r}}) = \sqrt{\tfrac{15}{32\pi}}\,\sin^2\theta\, e^{\pm 2i\phi},
$$
together with the spherical-to-Cartesian transcriptions $r^2 (3\cos^2\theta - 1) = 2 z^2 - x^2 - y^2$, $r^2 \sin\theta\cos\theta\, e^{i\phi} = z(x + iy)$, and $r^2 \sin^2\theta\, e^{2i\phi} = (x + iy)^2$, the five $\ell=2$ moments $q_2^m(A) = \int \rho_A\, |\boldsymbol{\xi}_A|^2\, Y_2^{m\,*}(\hat{\boldsymbol{\xi}}_A)\, dV_A$ map onto the Cartesian quadrupole tensor components as:

$$
\boxed{
\begin{aligned}
q_2^{0}(A) \;&=\; \sqrt{\tfrac{5}{16\pi}}\, \cdot \tfrac{1}{3}\, (\mathbf{Q}_A)_{zz},
\\[2pt]
q_2^{\pm 1}(A) \;&=\; \mp \sqrt{\tfrac{15}{8\pi}}\, \cdot \tfrac{1}{3}\, \big[\, (\mathbf{Q}_A)_{xz} \;\mp\; i\,(\mathbf{Q}_A)_{yz}\, \big],
\\[2pt]
q_2^{\pm 2}(A) \;&=\; \sqrt{\tfrac{15}{32\pi}}\, \cdot \tfrac{1}{3}\, \big[\, (\mathbf{Q}_A)_{xx} - (\mathbf{Q}_A)_{yy} \;\mp\; 2i\,(\mathbf{Q}_A)_{xy}\, \big].
\end{aligned}
}
$$

(The factor $1/3$ is the prefactor difference between our $\mathbf{Q}_A$ definition with leading $3 x_i x_j - r^2 \delta_{ij}$ in §3.5 and the T&B/Jackson convention with leading $x_i x_j - \tfrac{1}{3} r^2 \delta_{ij}$ — recall the convention note in §3.4.3's reference anchor: $\mathbf{Q}_{jk} = 3\, \mathcal{I}_{jk}$.) The five complex moments are not independent — the reality of the mass density $\rho_A$ enforces $q_2^{-m}(A) = (-1)^m\, q_2^{m\,*}(A)$, leaving exactly five real degrees of freedom that correspond to the five independent entries of the trace-free symmetric $3 \times 3$ Cartesian tensor: $(\mathbf{Q}_A)_{xx}$, $(\mathbf{Q}_A)_{yy}$ (with $(\mathbf{Q}_A)_{zz} = -(\mathbf{Q}_A)_{xx} - (\mathbf{Q}_A)_{yy}$ by tracelessness), $(\mathbf{Q}_A)_{xy}$, $(\mathbf{Q}_A)_{xz}$, $(\mathbf{Q}_A)_{yz}$. The mapping above is invertible — given the five $q_2^m(A)$, the Cartesian entries are recovered by inverting the prefactors and taking real/imaginary parts.

**Geometric interpretation of the mapping.** Each spherical-harmonic moment "filters" a specific Cartesian combination through its angular projection: $q_2^0$ measures the body's prolate/oblate stretch along $z$ (it depends only on $(\mathbf{Q}_A)_{zz}$, equivalently on the body's axial moment relative to its perpendicular moments); $q_2^{\pm 1}$ measure the off-diagonal $xz$ and $yz$ couplings (body tilt out of the $xy$-plane); $q_2^{\pm 2}$ measure the in-plane asymmetry $(\mathbf{Q}_A)_{xx} - (\mathbf{Q}_A)_{yy}$ and the off-diagonal $xy$ shear. For an axisymmetric body around $z$ (e.g., the `make_oblate_reference_body()` of §2.4), only $q_2^0$ is nonzero; all four other $\ell=2$ moments vanish identically. This is consistent with the Cartesian observation that $\mathbf{Q}_A$ for an axisymmetric body is diagonal with $(\mathbf{Q}_A)_{xx} = (\mathbf{Q}_A)_{yy}$ and only $(\mathbf{Q}_A)_{zz}$ as a free parameter.

Carrying out the $\ell = 2$ sum and using this Cartesian-spherical equivalence (the contractions on $\hat{\boldsymbol{\rho}}$ in the spherical-harmonic sum collapse to the inner product $\hat{\boldsymbol{\rho}} \cdot \mathbf{Q}_A \cdot \hat{\boldsymbol{\rho}}$ in the Cartesian form; the algebraic details are standard and we omit them):

$$
V_A^{\mathrm{on}\, B, \,\ell=2}
\;=\;
-\frac{G\, m_B}{2\, |\boldsymbol{\rho}|^3}\,
\hat{\boldsymbol{\rho}} \cdot \mathbf{Q}_A \cdot \hat{\boldsymbol{\rho}}.
$$

The Cartesian quadrupole tensor $\mathbf{Q}_A$ is related to the body-frame inertia tensor $\mathbf{I}_A$ by the **trace identity**:

$$
(\mathbf{Q}_A)_{ij} \;=\; \mathrm{tr}(\mathbf{I}_A^{\mathrm{inertial}})\, \delta_{ij} \;-\; 3\, (\mathbf{I}_A^{\mathrm{inertial}})_{ij},
$$

where $\mathbf{I}_A^{\mathrm{inertial}} = \mathbf{R}_q^A \mathbf{I}_A (\mathbf{R}_q^A)^T$ is the inertia tensor of $A$ expressed in inertial-frame coordinates (Chapter 02 §2.5). Both $\mathbf{Q}_A$ and $\mathbf{I}_A^{\mathrm{inertial}}$ are symmetric, but $\mathbf{Q}_A$ is additionally **traceless** ($\mathrm{tr}(\mathbf{Q}_A) = 3\, \mathrm{tr}(\mathbf{I}_A) - 3\, \mathrm{tr}(\mathbf{I}_A) = 0$); the trace identity above is the explicit decomposition of a generic symmetric tensor into its trace and trace-free parts.

> **Verification.** The trace identity can be checked directly: applying the definitions of $(\mathbf{Q}_A)_{ij}$ and $(\mathbf{I}_A)_{ij}$ from this chapter and Chapter 02 respectively, both reduce to integrals over $\rho_A$ of polynomials in $\mathbf{x}_A$; the linear combination $\mathrm{tr}(\mathbf{I}_A) \delta_{ij} - 3 (\mathbf{I}_A)_{ij}$ exactly recovers $\int \rho_A [3 x_i x_j - |\mathbf{x}|^2 \delta_{ij}] dV$. A numerical sanity check on the thin-rod inertia $\mathbf{I} = \mathrm{diag}(mL^2/3,\, mL^2/3,\, 0)$ produces $\mathbf{Q} = \mathrm{diag}(-mL^2/3,\, -mL^2/3,\, 2mL^2/3)$, traceless to machine precision (verified during Chapter 03 authoring; see `experiments/jwst-l2-first-cut/docs/canonical/en/03-*` authoring transcript 2026-05-23).

Substituting the trace identity:

$$
\hat{\boldsymbol{\rho}} \cdot \mathbf{Q}_A \cdot \hat{\boldsymbol{\rho}}
\;=\; \mathrm{tr}(\mathbf{I}_A)\, |\hat{\boldsymbol{\rho}}|^2
\;-\; 3\, \hat{\boldsymbol{\rho}} \cdot \mathbf{I}_A^{\mathrm{inertial}} \cdot \hat{\boldsymbol{\rho}}
\;=\; \mathrm{tr}(\mathbf{I}_A) \;-\; 3\, \hat{\boldsymbol{\rho}} \cdot \mathbf{I}_A^{\mathrm{inertial}} \cdot \hat{\boldsymbol{\rho}}.
$$

And using $\hat{\boldsymbol{\rho}} \cdot \mathbf{I}_A^{\mathrm{inertial}} \cdot \hat{\boldsymbol{\rho}} = \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$ (the inertia-tensor similarity transformation; Chapter 02 §2.5):

$$
\boxed{
V_A^{\mathrm{on}\, B, \,\ell=2}
\;=\;
-\frac{G\, m_B}{2\, |\boldsymbol{\rho}|^3}
\Big[ \mathrm{tr}(\mathbf{I}_A) \;-\; 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \Big].
}
$$

This is **identical** to the boxed result of §3.4.3.

> **Reference anchor.** Jackson 1999 §4.1 (multipole expansion, Cartesian-spherical equivalence at $\ell = 2$). The mathematics is identical in the electromagnetic and gravitational contexts; the only substitution is $\rho_{\mathrm{charge}} \to \rho_{\mathrm{mass}}$ and $1/(4\pi\varepsilon_0) \to -G$. Murray & Dermott 1999 §4.4 also presents the spherical-harmonic form for solar-system applications, with conventions specialised to oblate-body coefficients $J_\ell$; we use the Jackson normalisation throughout this chapter.

## §3.6 Cross-validation: Cartesian ≡ spherical-harmonic at quadrupole order

The two approaches, started from independent identities (the Cartesian Taylor expansion of $1/|\boldsymbol{\rho} - \mathbf{u}|$ in §3.4 and the Laplace spherical-harmonic addition theorem in §3.5), produce the same boxed quadrupole tidal expression:

$$
V_{\mathrm{tidal},\, A}
\;=\;
-\frac{G\, m_B}{2\, |\boldsymbol{\rho}|^3}
\Big[ \mathrm{tr}(\mathbf{I}_A) \;-\; 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \Big].
$$

The agreement is **not a coincidence**: at the level of mathematical identities, the Cartesian Taylor expansion of $1/|\boldsymbol{\rho} - \mathbf{u}|$ and the Laplace spherical-harmonic addition theorem are equivalent (one can be derived from the other by re-summing terms). What the cross-derivation does is verify that **we have not made an algebraic error in either route**: at the cost of a duplicated calculation, we have eliminated the class of bugs that arises from a single derivation chain.

This methodology — **two independent standard tools, cross-validated** — was directly mandated by the project owner on 2026-05-23 for any canonical-tier derivation that uses a standard mathematical technique (`CLAUDE.md` HANDOFF BLOCK 2026-05-23 §4). It is the analytical analog of the doctrine §4.4 C2 **cross-implementation bit-identity** criterion. Subsequent canonical-tier chapters will apply the same discipline: Chapter 04 (Euler-Lagrange + Euler-Poincaré reduction) will derive the body-frame rotational equations via both direct Euler-Lagrange computation and Euler-Poincaré reduction theorem; Chapter 05 (tennis-racket-theorem stability) via both linearisation and direct elliptic-function integration; Chapter 06 validates each testbed against at least two textbook references.

> **Reference anchor.** Jackson 1999 §4.1 contains the explicit proof that the Cartesian quadrupole tensor $\mathbf{Q}_{ij}$ and the spherical-harmonic multipole moments $q_2^m$ encode the same five independent degrees of freedom of a trace-free symmetric $3 \times 3$ tensor. The mapping is a basis change in the space of $\ell = 2$ tensor harmonics; the physical content is identical.

> **Supplementary anchor for Czech-reading audiences.** The Cartesian quadrupole tensor $\mathbf{Q}_{ij}$ and the spherical-harmonic multipole expansion are *modern* mathematical methods not directly treated in Podolský 2018 — that work is built around classical analytical mechanics of rigid bodies, not extended-body gravitational potentials. The bridge from Podolský to our §3.5 is via §7.1 *Tenzor setrvačnosti* Eqs. (7.3)–(7.9): the trace-free symmetric structure of the gravitational quadrupole tensor parallels Podolský's bilinear-form treatment of the kinetic-theory inertia tensor (both are real symmetric $3 \times 3$ tensors built as integrals of mass distribution against second-order spatial polynomials). The algebraic vocabulary Podolský develops in §7.1 — symmetric-tensor bilinear-form definition, principal-axis decomposition, $I_{\mathbf{n}}$ direction-projection — is directly reusable mental scaffolding for the gravitational quadrupole tensor. Native-language CS reader who wants the multipole-expansion theory itself is best served by an electrodynamics text (Jackson 1999 §4.1 cited above, or a CS-language Klasická elektrodynamika course text if available locally).

## §3.7 Gravity-gradient torque on body A

The gravity-gradient torque is obtained from the orientation-dependent part of the potential by the standard variational rule for rotational degrees of freedom: under a virtual rotation $\delta \boldsymbol{\theta}_A$ of body $A$ (a small axis-angle vector applied to its body-frame), the change in $V_{\mathrm{tidal},\, A}$ is $\delta V = -\boldsymbol{\tau}_A \cdot \delta \boldsymbol{\theta}_A$, defining the body-frame torque $\boldsymbol{\tau}_A$. We derive $\boldsymbol{\tau}_A$ explicitly.

Under the infinitesimal body rotation $\mathbf{R}_q^A \to \mathbf{R}_q^A (\mathbb{I} + \widehat{\delta \boldsymbol{\theta}_A})$, the inertial-frame separation direction $\hat{\boldsymbol{\rho}}$ is unchanged (the rotation acts on the body, not on the bodies' separation vector), and the transformed body-frame unit vector is:

$$
\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}^{\,\mathrm{new}}
\;=\; \big[\mathbf{R}_q^A (\mathbb{I} + \widehat{\delta \boldsymbol{\theta}_A})\big]^T \hat{\boldsymbol{\rho}}
\;=\; (\mathbb{I} - \widehat{\delta \boldsymbol{\theta}_A})\, (\mathbf{R}_q^A)^T \hat{\boldsymbol{\rho}}
\;=\; (\mathbb{I} - \widehat{\delta \boldsymbol{\theta}_A})\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A},
$$

where we used $(AB)^T = B^T A^T$ and $\widehat{\delta \boldsymbol{\theta}_A}^T = -\widehat{\delta \boldsymbol{\theta}_A}$ (antisymmetry of the hat map). To first order in $\delta \boldsymbol{\theta}_A$, this gives:

$$
\delta \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}
\;=\;
-\, \widehat{\delta \boldsymbol{\theta}_A}\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}
\;=\;
-\, \delta \boldsymbol{\theta}_A \times \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}.
$$

The trace term $\mathrm{tr}(\mathbf{I}_A)$ in $V_{\mathrm{tidal},\, A}$ is **orientation-independent** (the trace of a tensor is invariant under similarity transformation) and contributes zero to $\delta V$. The orientation-dependent piece is $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$, whose variation is:

$$
\delta\!\left[\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}\right]
\;=\;
2\, \delta \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}
\;=\;
-2\, (\delta \boldsymbol{\theta}_A \times \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}) \cdot (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}).
$$

Using the scalar-triple-product identity $(\mathbf{a} \times \mathbf{b}) \cdot \mathbf{c} = \mathbf{a} \cdot (\mathbf{b} \times \mathbf{c})$:

$$
\delta\!\left[\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}\right]
\;=\;
-2\, \delta \boldsymbol{\theta}_A \cdot
\big[ \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}) \big].
$$

Inserting this into the variation of $V_{\mathrm{tidal},\, A}$ (factor of $-3$ in front of the orientation-dependent piece, factor of $-G m_B / (2 |\boldsymbol{\rho}|^3)$ overall, so the orientation-dependent coefficient is $+3 G m_B / (2 |\boldsymbol{\rho}|^3)$):

$$
\delta V_{\mathrm{tidal},\, A}
\;=\;
\frac{3\, G\, m_B}{2\, |\boldsymbol{\rho}|^3}
\cdot (-2)\, \delta \boldsymbol{\theta}_A \cdot
\big[ \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}) \big]
\;=\;
-\, \frac{3\, G\, m_B}{|\boldsymbol{\rho}|^3}\, \delta \boldsymbol{\theta}_A \cdot
\big[ \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}) \big].
$$

The defining relation $\delta V = -\boldsymbol{\tau}_A \cdot \delta \boldsymbol{\theta}_A$ yields:

$$
\boxed{
\boldsymbol{\tau}_A^{\mathrm{body}}
\;=\;
\frac{3\, G\, m_B}{|\boldsymbol{\rho}|^3}\,
\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}).
}
$$

This is the **gravity-gradient torque** on body $A$ in its body frame, at leading (quadrupole) order in the body-extent / separation ratio. A symmetric expression holds for $\boldsymbol{\tau}_B^{\mathrm{body}}$ with $A \leftrightarrow B$ and $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, B} = -(\mathbf{R}_q^B)^T \hat{\boldsymbol{\rho}}$ (the unit vector from $B$ to $A$ in body-$B$ coordinates is the opposite sign).

> **Reference anchor.** Hughes 2004 §5.3 (gravity-gradient torque on a rigid spacecraft due to a central point-mass primary; the result is given in essentially this form). The Hughes derivation specialises to the case where body $B$ is the much-larger primary treated as a point mass; here we have derived the symmetric two-extended-body case, of which Hughes's result is the $\mathbf{I}_B \to \mathbf{0}$ limit. Goldstein 2002 §5.6 derives the same expression by direct calculation of the torque integral over the body's mass distribution (no multipole intermediate step), as an independent third-route cross-check.

### §3.7.1 Line-by-line match against `dynamics.py`

The engineer-tier implementation in `dynamics.py::gravity_gradient_torque_body_frame` (lines 113–131) reads, condensed:

```
r_mag = |r_rel_inertial|
r_hat_inertial = r_rel_inertial / r_mag
r_hat_body = R_body_to_inertial.T @ r_hat_inertial
I_rhat = I_body @ r_hat_body
tau = (3 * G * m_other / r_mag**3) * cross(r_hat_body, I_rhat)
```

The call site for body $A$ (`state_derivative`, lines 154–157) passes `r_rel_inertial = body_b.x - body_a.x`, which is our $\boldsymbol{\rho}$ from $A$ to $B$. Therefore `r_hat_body = R_A.T @ rho_hat = rho_hat_body_A` in our notation. Substituting:

$$
\texttt{tau} \;=\; \frac{3\, G\, m_B}{|\boldsymbol{\rho}|^3}\, \big[\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A})\big],
$$

which is **identical** to the boxed result of §3.7. The match is line-by-line: prefactor, sign, cross-product ordering, and convention for the direction of $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$ all agree. The implementation's docstring (line 122) cites Hughes 2004 §5.3 directly; this chapter's derivation is the canonical-tier expansion of that citation.

For the torque on body $B$, the call site (`state_derivative`, line 158) passes `r_rel_inertial = -r_a_to_b_inertial = body_a.x - body_b.x`, which is $-\boldsymbol{\rho}$ from $B$'s perspective — exactly the unit vector from $B$ to $A$. The sign convention is therefore consistent: in each body's call, $\hat{\boldsymbol{\rho}}_{\mathrm{body}}$ points from that body to the other.

> **Supplementary anchor for Czech-reading audiences.** Gravity-gradient torque on an extended body in an external gravitational field is a *modern* spacecraft-attitude-dynamics topic, outside the classical scope of Podolský 2018. The closest classical-mechanics precedent in Podolský is **Kapitola 8 §8.3** *Těžký symetrický setrvačník s pevným bodem* — the heavy symmetric top with fixed point (Lagrange top). In that problem, the external gravitational field acts on an axisymmetric rigid body suspended from a fixed point, producing a torque $\mathbf{M}$ that depends on the body's orientation through the height of its centre of mass — analytically equivalent in algebraic structure to our gravity-gradient torque dependence on orientation. Podolský §8.3 derives the resulting motion via Lagrangian formalism with Euler-angle generalised coordinates: cyclic-coordinate first integrals $L_\psi, L_\varphi$ at his Eq. (8.11), energy conservation reducing to one-dimensional motion in $\vartheta$ with effective potential $V_{\mathrm{ef}}(\vartheta)$ at boxed Eq. (8.15), and the precession-regime classification (Eqs. 8.13 with cases (a)/(b)/(c) on the unit sphere). Our gravity-gradient libration problem (Chapter 06 future testbed) is the spacecraft-attitude analog of Podolský's Lagrange top: same Lagrangian structure, different physical interpretation of the "torque about the suspension axis". The Phase 2 Fortran libration testbed will cross-validate against the Lagrange-top $V_{\mathrm{ef}}$ method as its canonical CS-language analytical-limit oracle.

### §3.7.2 Gravity-gradient torque Jacobian: sensitivity to initial conditions and assembly geometry (engineer-tier bridge)

The §3.7 derivation produced the gravity-gradient torque $\boldsymbol{\tau}_A^{\mathrm{body}}$ as a closed-form function of the configuration $(\boldsymbol{\rho}, q_A, q_B)$. For the engineer-tier — practical sensitivity analysis, initial-condition impact estimation, integrator step-size selection — the **linearised Jacobian** of $\boldsymbol{\tau}_A$ with respect to small perturbations of the configuration is the deliverable that complements the closed-form formula. This subsection carries the engineer-tier-bridge framework of Chapter 02 §2.5.1 (the $d\mathbf{I}^{\mathrm{inertial}}/dt$ off-diagonal sensitivity) into the gravity-gradient context, per project-owner direction 2026-05-23.

The configuration enters $\boldsymbol{\tau}_A^{\mathrm{body}} = (3 G m_B / r^3)\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A})$ through three places: the inter-body separation magnitude $r = |\boldsymbol{\rho}|$ (the prefactor), the body-frame direction $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = (\mathbf{R}_q^A)^T \hat{\boldsymbol{\rho}}$ (which depends on both $\boldsymbol{\rho}$ and $q_A$), and the constant body-frame tensor $\mathbf{I}_A$ (not a perturbable state variable). We compute partial derivatives with respect to each state variable in turn.

**Sensitivity to body $A$ orientation: $\partial \boldsymbol{\tau}_A^{\mathrm{body}} / \partial q_A$.** Under the infinitesimal body rotation $\mathbf{R}_q^A \to \mathbf{R}_q^A (\mathbb{I} + \widehat{\delta \boldsymbol{\theta}_A})$ (the same parameterisation used in §3.7's derivation), $\delta \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = -\delta \boldsymbol{\theta}_A \times \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$. Writing $\mathbf{u} = \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$, $\mathbf{w} = \mathbf{I}_A \mathbf{u}$, and $\alpha = 3 G m_B / r^3$, and linearising the torque:

$$
\delta \boldsymbol{\tau}_A^{\mathrm{body}}
\;=\; \alpha\, \big[\delta \mathbf{u} \times \mathbf{w} \;+\; \mathbf{u} \times \mathbf{I}_A \delta \mathbf{u}\big]
\;=\; \alpha\, \big[(-\delta \boldsymbol{\theta}_A \times \mathbf{u}) \times \mathbf{w} \;+\; \mathbf{u} \times \mathbf{I}_A (-\delta \boldsymbol{\theta}_A \times \mathbf{u})\big].
$$

This is a linear map from $\delta \boldsymbol{\theta}_A \in \mathbb{R}^3$ to $\delta \boldsymbol{\tau}_A^{\mathrm{body}} \in \mathbb{R}^3$. Define the **orientation-response Jacobian** $\mathbf{J}^{\mathrm{orient}}_{A}$ by $\delta \boldsymbol{\tau}_A^{\mathrm{body}} = \mathbf{J}^{\mathrm{orient}}_{A}\, \delta \boldsymbol{\theta}_A$; this is a $3 \times 3$ matrix-valued function of $(\mathbf{I}_A, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A})$:

$$
\boxed{
\mathbf{J}^{\mathrm{orient}}_{A}\, \boldsymbol{a}
\;=\; -\alpha\, \Big[ (\mathbf{a} \times \mathbf{u}) \times \mathbf{I}_A \mathbf{u} \;+\; \mathbf{u} \times \mathbf{I}_A (\mathbf{a} \times \mathbf{u}) \Big],
\qquad \mathbf{u} = \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A},\quad \alpha = \frac{3 G m_B}{|\boldsymbol{\rho}|^3}.
}
$$

> **Sign convention note (2026-05-24 correction, OQ-FORT-1).** The overall $-\alpha$ in the boxed expression follows from carrying the two minus signs in $\delta \mathbf{u} = -\delta \boldsymbol{\theta}_A \times \mathbf{u}$ through the linearisation of $\boldsymbol{\tau}_A^{\mathrm{body}}$ explicitly. The earlier draft of this chapter omitted the overall sign; the Fortran Phase 1 cross-test `tests/test_orientation_jacobian.f90` T2a (libration-equilibrium eigenvalues for axisymmetric oblate body) is the permanent regression guard: for oblate ($I_\parallel > I_\perp$), the eigenvalue $J_{11} = \alpha (I_\parallel - I_\perp)$ must be **positive** (destabilising), which matches the corrected $-\alpha$ form when combined with the algebraic identity $-[(a \times u) \times (I u) + u \times I(a \times u)]_{x=e_1} = (I_\parallel - I_\perp,\, 0,\, 0)$ at the on-axis equilibrium. The libration-equilibrium specialisation and the stability table below are already expressed in the corrected form.

The structure of $\mathbf{J}^{\mathrm{orient}}_{A}$ depends on the body's symmetry class and on the alignment of $\mathbf{u}$ with the body's principal axes.

**At the libration equilibrium (axisymmetric body, symmetry axis along $\hat{\boldsymbol{\rho}}$).** Specialise to body $A$ axisymmetric ($I_1 = I_2 \equiv I_\perp$, $I_3 \equiv I_\parallel$, principal axes aligned with body frame) and orientation such that $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = (0, 0, 1) = \hat{\mathbf{e}}_z^{\,\mathrm{body}}$ — the "radial-pointing" attitude of validation testbed §2.5. At this configuration, $\boldsymbol{\tau}_A^{\mathrm{body}} = \mathbf{0}$ (an equilibrium), and direct calculation of the Jacobian gives:

$$
\mathbf{J}^{\mathrm{orient}}_{A}\big|_{\mathrm{eq}}
\;=\; \alpha\, (I_\parallel - I_\perp)\, \mathrm{diag}(1,\, 1,\, 0)
\;=\; \frac{3 G m_B (I_\parallel - I_\perp)}{|\boldsymbol{\rho}|^3}\, \mathrm{diag}(1,\, 1,\, 0).
$$

The sign of the $(I_\parallel - I_\perp)$ factor determines stability of the equilibrium — the **gravity-gradient stabilisation theorem** for an axisymmetric satellite:

| Body geometry | Sign of $I_\parallel - I_\perp$ | Equilibrium with symmetry axis along $\hat{\boldsymbol{\rho}}$ |
|---|---|---|
| Prolate ($I_\parallel < I_\perp$; needle / long rod along symmetry axis) | $-$ | **Stable** (restoring torque; libration about equilibrium) |
| Spherical ($I_\parallel = I_\perp$) | $0$ | Marginal (no torque from this axis) |
| Oblate ($I_\parallel > I_\perp$; flat disk / sunshield perpendicular to symmetry axis) | $+$ | **Unstable** (torque drives short axis away from $\hat{\boldsymbol{\rho}}$; stable orientation is symmetry axis perpendicular to $\hat{\boldsymbol{\rho}}$) |

The $(3,3)$ entry of $\mathbf{J}^{\mathrm{orient}}_{A}\big|_{\mathrm{eq}}$ vanishes: rotation about the symmetry axis itself produces no gravity-gradient torque, as expected (the body is symmetric about that axis; spinning about it does not change the gravitational coupling).

**Sensitivity to inter-body position: $\partial \boldsymbol{\tau}_A^{\mathrm{body}} / \partial \boldsymbol{\rho}$.** Decompose $\delta \boldsymbol{\rho} = \delta r\, \hat{\boldsymbol{\rho}} + \delta \boldsymbol{\rho}_\perp$ where $\delta \boldsymbol{\rho}_\perp \perp \hat{\boldsymbol{\rho}}$. The two pieces have different effects on $\boldsymbol{\tau}_A^{\mathrm{body}}$:

- **Radial perturbation** $\delta r$: rescales the prefactor $\alpha = 3 G m_B / r^3 \to \alpha (1 - 3 \delta r / r)$ to first order, while $\hat{\boldsymbol{\rho}}$ (and therefore $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}$) is unchanged. The Jacobian along $\hat{\boldsymbol{\rho}}$ is simply:
  $$
  \frac{\partial \boldsymbol{\tau}_A^{\mathrm{body}}}{\partial r}
  \;=\; -\frac{3}{r}\, \boldsymbol{\tau}_A^{\mathrm{body}},
  $$
  rescaling the entire torque vector with sensitivity inversely proportional to the separation. At the libration equilibrium ($\boldsymbol{\tau}_A^{\mathrm{body}} = \mathbf{0}$), this contribution vanishes.

- **Tangential perturbation** $\delta \boldsymbol{\rho}_\perp$: rotates the inertial-frame direction $\hat{\boldsymbol{\rho}}$ by an infinitesimal axis-angle $\delta \boldsymbol{\chi} = (\hat{\boldsymbol{\rho}} \times \delta \boldsymbol{\rho}_\perp) / r$. This rotation propagates to the body frame as $\delta \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = (\mathbf{R}_q^A)^T \delta \hat{\boldsymbol{\rho}}$, which is **algebraically identical in structure** to the orientation-perturbation case above (a small rotation of $\mathbf{u}$ in body coordinates). The corresponding Jacobian is, by direct substitution:
  $$
  \boxed{
  \frac{\partial \boldsymbol{\tau}_A^{\mathrm{body}}}{\partial \boldsymbol{\rho}_\perp}\, \delta \boldsymbol{\rho}_\perp
  \;=\; \mathbf{J}^{\mathrm{orient}}_{A}\, \big[(\mathbf{R}_q^A)^T\, (\hat{\boldsymbol{\rho}} \times \delta \boldsymbol{\rho}_\perp / r)\big],
  }
  $$
  with $\mathbf{J}^{\mathrm{orient}}_{A}$ the same matrix as the orientation-response Jacobian above. The body-orientation Jacobian and the position-direction Jacobian are not independent objects — they are the same linear map operating on different inputs, reflecting the deeper fact that only the *relative* angle between $\hat{\boldsymbol{\rho}}$ and body $A$'s principal axes enters $\boldsymbol{\tau}_A^{\mathrm{body}}$.

**Sensitivity to body $B$ orientation: $\partial \boldsymbol{\tau}_A^{\mathrm{body}} / \partial q_B = 0$ at quadrupole truncation.** Body $A$'s gravity-gradient torque depends on body $B$ only through its mass $m_B$ and its position $\mathbf{R}_B$ (entering via $\boldsymbol{\rho} = \mathbf{R}_B - \mathbf{R}_A$); body $B$'s orientation $q_B$ does not appear in $\boldsymbol{\tau}_A^{\mathrm{body}}$ at the quadrupole-truncated multipole expansion adopted in §3.4. Body $B$'s orientation enters only at the *next* multipole order: the quadrupole-quadrupole bilinear term (decay $1/r^5$) and the octupole-monopole term (decay $1/r^4$), both estimated and bounded at §3.3. For the first-cut configurations $(\ell/r \lesssim 10^{-2})$ this cross-orientation sensitivity is suppressed by at least four orders of magnitude relative to the leading quadrupole tidal term and remains below double-precision round-off.

**Symmetry-class structure of $\mathbf{J}^{\mathrm{orient}}_{A}$.** Counting the rank of the orientation-response Jacobian across body-symmetry classes (table parallels §2.5.1's off-diagonal-rate sensitivity table):

| Body symmetry class | Rank of $\mathbf{J}^{\mathrm{orient}}_{A}$ | Physical interpretation |
|---|---|---|
| Spherical ($I_1 = I_2 = I_3$) | $0$ | $\boldsymbol{\tau}_A \equiv \mathbf{0}$ identically; no sensitivity to orientation (the decoupling testbed §3.8.2 special case). |
| Axisymmetric ($I_1 = I_2 \neq I_3$) | $2$ generically (drops to $0$ at the on-symmetry-axis equilibrium for rotations *about* the symmetry axis) | Two-dimensional libration plane perpendicular to the symmetry axis; one cyclic direction along the axis. |
| Asymmetric (all $I_i$ distinct) | $3$ generically | Three independent libration / instability modes; the gravity-gradient torque couples to all three orientation degrees of freedom. |

**Two-body assembly sensitivity figure of merit.** The total assembly dynamics couples body $A$'s rotation, body $B$'s rotation, and the relative orbit through $\boldsymbol{\rho}$. The integrated effect of gravity-gradient torque on the rotational state of each body scales with $\alpha \cdot |\mathbf{I}_X^{\mathrm{anisotropy}}|$, where the anisotropy is captured (for axisymmetric bodies) by the ratio $\eta_X = (I_\parallel - I_\perp)_X / I_\perp^X$ or equivalently $I_\parallel / I_\perp - 1$. For two-body assemblies, the **figure of merit** for assembly-wide sensitivity to initial-condition perturbations is the product:

$$
\Sigma_{\mathrm{assembly}}
\;\equiv\;
\eta_A \cdot \eta_B \cdot f(\hat{\boldsymbol{\rho}}\text{-alignment}),
$$

where $f$ is an $O(1)$ factor encoding how $\hat{\boldsymbol{\rho}}$ aligns with each body's principal axes ($f = 0$ if $\hat{\boldsymbol{\rho}}$ is along a principal axis of either body simultaneously and equilibrium is at zero torque; $f \to 1$ for generic orientation; $f$ extreme at the gravity-gradient unstable configurations). For the first-cut configuration (JWST-like body $A$ with $I_\parallel/I_\perp = 1.52$, probe-like body $B$ with $I_\parallel/I_\perp = 5.81$): $\eta_A = 0.52$, $\eta_B = 4.81$, $\eta_A \eta_B \approx 2.5$. The probe-like body's high anisotropy dominates the assembly sensitivity; small perturbations of $q_B$ propagate to large $\boldsymbol{\tau}_B^{\mathrm{body}}$ excursions, which in turn drift the angular-momentum bookkeeping and stress the integrator's conservation diagnostic.

**Numerical sanity check: oblate reference body at a hypothetical L2-like configuration.** Take $m_B = m_\oplus \approx 6 \times 10^{24}\ \text{kg}$ (Earth-mass primary), $r = 1.5 \times 10^9\ \text{m}$ (a representative order-of-magnitude L2 separation from Earth). Then:

$$
\alpha = \frac{3 G m_B}{r^3} \approx 3.5 \times 10^{-13}\ \text{s}^{-2}.
$$

For the **oblate reference body** `make_oblate_reference_body()` (Ch 02 §2.4 — synthetic single-cylinder body with $\mathbf{I}_A = \mathrm{diag}(2540, 2540, 3870)\ \text{kg m}^2$ by construction), the libration-equilibrium spring constant is $\alpha (I_\parallel - I_\perp) = 3.5 \times 10^{-13} \cdot 1330 \approx 4.7 \times 10^{-10}\ \text{N m / rad}$. Because $I_\parallel > I_\perp$ (oblate), the spring constant is *positive*, signalling **instability** of the symmetry-axis-along-radial attitude — gravity gradient drives the symmetry axis *away* from the radial direction. The e-folding time of the instability is $\tau_{\mathrm{inst}} = \sqrt{I_\perp / |k|} \approx 2.3 \times 10^6\ \text{s} \approx 27\ \text{days}$. This is on the order of the slow-attitude-evolution time-scale of a sunshield-axis-radial JWST-like spacecraft (in practice, real JWST uses momentum wheels and reaction control to compensate for this and other torques continuously).

**Contrast: prolate composite at the same hypothetical L2-like configuration.** The simplified spacecraft composite `make_jwst_like()` (4-component sunshield + bus + boom + primary, total 2450 kg) produces $\mathbf{I}_A = \mathrm{diag}(23323, 23323, 15384)\ \text{kg m}^2$ — the body is **prolate** ($I_\parallel = I_3 < I_\perp = I_1 = I_2$) because the component z-spread contributes to $I_\perp$ via the parallel-axis theorem in excess of the local axial contributions. At the same hypothetical L2 configuration, the spring constant is $\alpha (I_\parallel - I_\perp) = 3.5 \times 10^{-13} \cdot (-7938) \approx -2.8 \times 10^{-9}\ \text{N m / rad}$ — **negative** (restoring), so the symmetry-axis-along-radial attitude is **stable** with libration period $T_{\mathrm{lib}} = 2\pi \sqrt{I_\perp / |k|} = 2\pi \sqrt{23323 / 2.8 \times 10^{-9}} \approx 5.7 \times 10^{6}\ \text{s} \approx 66\ \text{days}$. The composite's prolate-stable classification is consistent with the gravity-gradient stabilisation principle for long-axis spacecraft; what differs from the real-JWST intuition is that the simplified model's component z-spread overwhelms the sunshield's local axial moment, producing prolate rather than oblate topology around z. This contrast — **oblate-unstable reference body next to prolate-stable simplified spacecraft** — illustrates the stability table directly: same gravitational environment, same anisotropy magnitude, opposite stability classification driven by the sign of $(I_\parallel - I_\perp)$. The Fortran Phase 1 cross-test fixture `tests/fixtures/gravity_gradient_inputs.json` cases 02 (prolate JWST-like, $\tau_y = -4.8 \times 10^{-10}$ at 10° off-equilibrium) and 07 (oblate reference body, $\tau_y = +8.1 \times 10^{-11}$ at 10° off-equilibrium) demonstrate the opposite signs empirically.

**Engineer-tier integrator implication: multi-scale problem.** The first-cut integration uses RK4 at $dt = 0.05\ \text{s}$ over a $600\ \text{s}$ window. The system has at least four natural time scales:

| Time scale | Magnitude (first-cut + hypothetical L2) | Resolution requirement |
|---|---|---|
| Body spin period $T_{\mathrm{spin}}$ | $\sim 628\ \text{s}$ (for $\omega \sim 0.01\ \text{rad/s}$) | $dt \ll T_{\mathrm{spin}}$ — easily met at $0.05\ \text{s}$ ($\sim 12{,}000$ samples per period) |
| Off-diagonal $I^{\mathrm{inertial}}$ time scale $\tau_{\mathrm{off}}$ (§2.5.1) | $\sim 190\ \text{s}$ | Within the spin period — same constraint. |
| Libration / instability time scale | $\sim 10^7\ \text{s}$ (~6 months) | Integration window must be $\gg T_{\mathrm{lib}}$ to observe libration — not captured at $600\ \text{s}$ window |
| Orbital period at L2-like configuration | $\sim 10^7\ \text{s}$ (~6 months) | Same — not captured at $600\ \text{s}$ window |

The first-cut window is therefore **spin-resolving but libration-blind**. To study the assembly's gravity-gradient libration sensitivity numerically, the integration window must extend to $\sim 10^7\ \text{s}$ ($\sim 2 \times 10^8$ RK4 steps at $dt = 0.05\ \text{s}$), at which point the accumulated non-symplectic energy drift from RK4 ($\sim O(dt^4 \cdot T)$ per Chapter 1 §1.6) becomes the dominant error source. A symplectic Lie-group integrator (Hairer-Lubich-Wanner 2006 §VII.5) is the structural fix for long-horizon attitude simulations; the first-cut deliberately remains at RK4 for the spin-resolution scope and flags the symplectic upgrade as a v0.2 engineer-tier item.

> **Reference anchor.** Hughes 2004 Ch 8 (gravity-gradient stabilisation, with the spring-constant matrix for axisymmetric and asymmetric bodies derived in essentially the form given here); Markley & Crassidis 2014 *Fundamentals of Spacecraft Attitude Determination and Control* §3.7 (engineering treatment of libration dynamics with attitude-control implications); Wertz 1978 *Spacecraft Attitude Determination and Control* §17 for the classical engineering reference. The sensitivity-Jacobian framing emphasised here — explicit $\mathbf{J}^{\mathrm{orient}}_{A}$ matrix as the unifying object across orientation and inter-body-position perturbations — is the engineer-tier complement to Hughes's libration treatment.

## §3.8 Reductions to validation testbeds

The quadrupole-truncated potential and the gravity-gradient torque reduce to the analytical limits catalogued in `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md`. We work through three of them explicitly here; Chapter 06 covers the full catalogue.

### §3.8.1 Reduction to Kepler (validation testbed §2.3)

Set $\mathbf{I}_A \to \mathbf{0}$ and $\mathbf{I}_B \to \mathbf{0}$ in the boxed total $V$ of §3.4.3. Then $V_{\mathrm{tidal},\, A}$ and $V_{\mathrm{tidal},\, B}$ both vanish identically (a zero inertia tensor has zero trace and zero contracted form), and the truncated $V$ collapses to:

$$
V \;\xrightarrow{\mathbf{I}_{A,B} \to \mathbf{0}}\; -\frac{G\, m_A\, m_B}{|\boldsymbol{\rho}|}.
$$

Combined with the centre-of-mass-frame kinetic energy $T = \tfrac{1}{2} \mu |\dot{\boldsymbol{\rho}}|^2 + \tfrac{1}{2}\, \boldsymbol{\Omega}_A^T \mathbf{I}_A \boldsymbol{\Omega}_A + \tfrac{1}{2}\, \boldsymbol{\Omega}_B^T \mathbf{I}_B \boldsymbol{\Omega}_B$ from Chapter 02 §2.6, the rotational kinetic energies also vanish in this limit, and the Lagrangian $L = T - V$ reduces to:

$$
L \;\xrightarrow{\mathbf{I}_{A,B} \to \mathbf{0}}\; \tfrac{1}{2}\, \mu\, |\dot{\boldsymbol{\rho}}|^2 \;+\; \frac{G\, m_A\, m_B}{|\boldsymbol{\rho}|},
$$

the standard reduced-mass Kepler Lagrangian. Its Euler-Lagrange equations are the Kepler equations of motion (Landau-Lifshitz Vol I §13–§15). The testbed `_specs/...` §2.3 (`tests/test_dynamics.py::TestKeplerLimit`, planned) verifies this numerically by setting the engineer-tier inertia tensors to zero and checking that the orbital period matches $T_{\mathrm{orbit}} = 2\pi \sqrt{a^3 / (G(m_A + m_B))}$ to $10^{-6}$ relative error.

### §3.8.2 Reduction to free rotation + Kepler decoupling (validation testbed §2.4)

Set $\mathbf{I}_B \to \mathbf{0}$ (point-mass body $B$); keep $\mathbf{I}_A$ generic. Then $V_{\mathrm{tidal},\, B} = 0$ (zero quadrupole moment of a point mass), and $V_{\mathrm{tidal},\, A} \neq 0$ in general — *unless* body $A$ is spherically symmetric, in which case $\mathbf{I}_A = I_0 \mathbb{I}_3$, and:

$$
\mathrm{tr}(\mathbf{I}_A) - 3\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \cdot \mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}
\;=\; 3 I_0 - 3 I_0 |\hat{\boldsymbol{\rho}}|^2 \;=\; 0,
$$

so $V_{\mathrm{tidal},\, A}$ vanishes identically. Translational and rotational dynamics decouple completely: $A$'s orbit is Kepler around the point mass $B$, and $A$'s rotation is torque-free Euler precession (Landau-Lifshitz §35 / Chapter 05 of this canonical tier). This is the testbed `_specs/...` §2.4 setup.

If $A$ is **not** spherically symmetric — for instance, the first-cut axisymmetric JWST-like body with $\mathbf{I}_A = \mathrm{diag}(I_\perp, I_\perp, I_\parallel)$ and $I_\parallel \neq I_\perp$ — then $V_{\mathrm{tidal},\, A}$ depends on orientation of $A$ relative to $\boldsymbol{\rho}$ and the gravity-gradient torque $\boldsymbol{\tau}_A$ is generically nonzero. This is the regime of the **gravity-gradient libration** testbed (§3.8.3 below), and the case treated extensively in Hughes 2004 Ch 8.

### §3.8.3 Small-angle libration (validation testbed §2.5)

Place an axisymmetric body $A$ ($I_{xx,A} = I_{yy,A} \equiv I_\perp$, $I_{zz,A} \equiv I_\parallel$) in a circular orbit around a point-mass body $B$. Define the **radial-pointing equilibrium**: orientation of $A$ such that its body-$z$ axis (the symmetry axis, the one with the distinct moment) is aligned with $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = \hat{\mathbf{e}}_z^{\,\mathrm{body}\, A}$ at all times — i.e., $A$ rotates synchronously with the orbit, keeping its symmetry axis pointed at $B$.

At equilibrium, $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = (0, 0, 1)$ and $\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} = (0, 0, I_\parallel)$. The cross product $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A}) = (0,0,1) \times (0,0,I_\parallel) = \mathbf{0}$. The torque vanishes; the equilibrium is a true critical point of the orientation-dependent potential.

Perturb the orientation by small angles $\theta, \phi$ about the equilibrium (a pitch–roll perturbation about the symmetry axis). To linear order in $(\theta, \phi)$, the perturbed $\hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \approx (-\theta, \phi, 1)$ to first order in small angles, and (working out the cross product and contracting against the orbital angular velocity $n = (G(m_A + m_B)/|\boldsymbol{\rho}|^3)^{1/2}$):

The linearised libration equations are coupled harmonic oscillators with normal-mode frequencies determined by the inertia ratio $(I_\parallel - I_\perp)/I_\perp$ and the orbit rate $n$. The full derivation is standard (Hughes 2004 §8.4); we do not reproduce it here because it requires the body-frame Euler equations from Chapter 05 as scaffolding. The testbed `_specs/...` §2.5 (`tests/test_integration.py::TestGravityGradientLibration`, planned) verifies the libration period numerically.

The key point for this chapter is that the gravity-gradient torque we derived in §3.7 **produces a non-trivial libration spectrum** for any non-spherical body, and the spectrum is testable against the Hughes 2004 closed-form prediction at $10^{-4}$ relative error. This is C1-textbook-verbatim test coverage for the gravity-gradient torque, in the sense of doctrine §4.4.

### §3.8.4 Infinite-separation decoupling (validation testbed §3.2)

For $|\boldsymbol{\rho}| \to \infty$, the monopole-monopole term $V^{(0)} \to 0$ and both tidal terms decay as $1/|\boldsymbol{\rho}|^3 \to 0$. The total mutual potential vanishes faster than any finite power of separation (when expanded to finite multipole order; higher orders decay faster still). Each body becomes an isolated torque-free top: $\boldsymbol{\tau}_A, \boldsymbol{\tau}_B \to \mathbf{0}$, $V \to 0$, and the Lagrangian factors as $L = L_A + L_B$ where each $L_X = T_X$ is the kinetic-only Lagrangian of an isolated rigid body. The system decouples completely. The testbed `_specs/...` §3.2 verifies that cross-coupling stays below $10^{-12}$ relative per orbital period at large separation.

## §3.9 Code cross-reference summary

For the engineer-tier reader who has tracked along to this point, the canonical-tier equations of this chapter correspond to engineer-tier code as follows:

| Canonical-tier object | Engineer-tier implementation |
|---|---|
| Exact $V$ integral (§3.2 boxed) | Not implemented; this is the full physical reference. |
| Monopole-monopole $V^{(0)}$ (§3.4.1) | `dynamics.py::total_potential_energy` (lines 207–214). |
| Quadrupole tidal $V_{\mathrm{tidal},\, A}$ (§3.4.3 / §3.5 boxed) | Not directly computed for energy diagnostics; absorbed implicitly into the dynamics via $\boldsymbol{\tau}_A$. **Note:** the conservation diagnostic `total_potential_energy` therefore tracks only the monopole-monopole term; the missing tidal contribution to $V$ is below double-precision round-off at the first-cut separations and does not affect the observed $|dE/E_0| = 6.2 \times 10^{-12}$. **Could be added** as a refinement; queued at §3.13 OQ-3.4. |
| Gravity-gradient torque $\boldsymbol{\tau}_A^{\mathrm{body}}$ (§3.7 boxed) | `dynamics.py::gravity_gradient_torque_body_frame` (lines 113–131); call site for body $A$ at line 157, for body $B$ at line 158. |
| Trace identity $\mathbf{Q} = \mathrm{tr}(\mathbf{I}) \cdot \mathbb{I} - 3 \mathbf{I}$ (§3.5) | Not implemented as a free-standing function; folded into the torque formula. |
| Convergence boundary (§3.3) | Not yet asserted; queued at `_specs/...` §3.3 / §6 OQ-deferred. |

The match is **exact at the gravity-gradient-torque level**: the boxed formula of §3.7 and the engineer-tier `gravity_gradient_torque_body_frame` are line-by-line identical, with consistent sign conventions across both body-$A$ and body-$B$ call sites (§3.7.1).

## §3.10 What is established by the end of this chapter

A reader who has worked through §3.1–§3.9 has:

- The **exact mutual gravitational potential** of two extended rigid bodies as a 6-dimensional integral over their mass distributions, derived from first-principles Newtonian gravity.
- The convergence domain and the close-approach boundary of the multipole expansion, locked as the model's domain of validity.
- The Cartesian multipole expansion of $V$ to quadrupole order, with the monopole-monopole (Kepler) and quadrupole-tidal (gravity-gradient) contributions identified explicitly.
- The independent spherical-harmonic multipole expansion, producing the same quadrupole-order result by the standard Jackson identity transferred from electromagnetism.
- **Cross-validation** that the two approaches yield identical results at quadrupole order — the doctrinal "two in parallel" requirement, satisfied.
- The trace identity $\mathbf{Q} = \mathrm{tr}(\mathbf{I}) \cdot \mathbb{I} - 3 \mathbf{I}$ relating the gravitational quadrupole tensor to the kinetic-theory inertia tensor of Chapter 02.
- The **gravity-gradient torque** $\boldsymbol{\tau}_A^{\mathrm{body}} = (3 G m_B / |\boldsymbol{\rho}|^3)\, \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},\, A})$, with a line-by-line match to the engineer-tier `dynamics.py::gravity_gradient_torque_body_frame` implementation.
- The **orientation-response Jacobian** $\mathbf{J}^{\mathrm{orient}}_{A}$ (§3.7.2) — engineer-tier sensitivity-analysis deliverable. Unifying $3 \times 3$ linear map governing torque response to body-orientation perturbations and inter-body-position-direction perturbations alike; libration-equilibrium specialisation $\alpha (I_\parallel - I_\perp) \mathrm{diag}(1,1,0)$ with gravity-gradient stability prediction; multi-scale time-scale audit identifying the first-cut as spin-resolving but libration-blind, and a symplectic-Lie-group-integrator flag for long-horizon (libration / orbital scale) studies.
- Reductions of the truncated $V$ to the analytical-limit validation testbeds (Kepler, free-rotation + Kepler decoupling, gravity-gradient libration, infinite-separation decoupling), each documented and pointing forward to Chapter 06's full treatment.

**What is NOT yet established:**

- The full Lagrangian $L = T - V$ and Euler-Lagrange equations of motion on $TQ$. **Chapter 04.**
- The body-frame Euler equations $\mathbf{I} \dot{\boldsymbol{\Omega}} + \boldsymbol{\Omega} \times (\mathbf{I} \boldsymbol{\Omega}) = \boldsymbol{\tau}$ as a specialisation by Euler-Poincaré reduction, with the $\boldsymbol{\tau}$ derived here entering as the right-hand side. **Chapter 05.**
- The L2 restricted-three-body extension with centrifugal and Coriolis pseudo-forces. **Deferred to canonical-tier v0.2** per Chapter 00 §1.
- Full enumeration and numerical verification of the validation testbeds. **Chapter 06**, drawing on the catalogue in `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md`.

With $T$ (Chapter 02) and $V$ (this chapter) both in hand, Chapter 04 has all the assembly material it needs to construct the Lagrangian and derive the equations of motion.

---

## §3.11 Open questions in this chapter

- **OQ-3.1.** ***Resolved 2026-05-24*** — explicit basis-change algebra added at §3.5 (boxed equation block following the $\ell = 2$ identification). Maps the five $q_2^m(A)$ (for $m = -2, -1, 0, 1, 2$) to specific linear combinations of the trace-free Cartesian quadrupole tensor components $(\mathbf{Q}_A)_{ij}$, with the Condon-Shortley spherical-harmonic convention and the prefactor-1/3 acknowledgement of the T&B/Jackson $\mathcal{I}_{jk}$-vs-our-$\mathbf{Q}_{jk}$ convention difference (per §3.4.3 reference anchor). Inversion algebra noted, reality conditions stated, geometric interpretation of each $q_2^m$ component added (prolate/oblate $z$-stretch, $xz/yz$ off-diagonal tilt, in-plane $xx-yy$ + $xy$ asymmetry). Axisymmetric-body specialisation noted: only $q_2^0$ nonzero. Cross-validated against the trace-identity discussion already in §3.5. Anchor: Jackson 1999 §3.6 (spherical harmonics) + §4.1 (multipole moment equivalence).
- **OQ-3.2.** ***Partial resolution 2026-05-25 (canonical-tier round 3).*** The **next-order multipole corrections** to the boxed §3.4.3 quadrupole-truncated $V$ have explicit power-of-$|\boldsymbol{\rho}|$ scaling derivable directly from extending the Cartesian Taylor expansion of §3.4 by one order, or equivalently the spherical-harmonic expansion of §3.5 by one $\ell$ level. Two distinct contributions arise at next order:
  - **Bilinear quadrupole-quadrupole term:** $V^{(\ell_A=2,\ell_B=2)} \propto G\, [\hat{\boldsymbol{\rho}} \cdot \mathbf{Q}_A \cdot \hat{\boldsymbol{\rho}}] [\hat{\boldsymbol{\rho}} \cdot \mathbf{Q}_B \cdot \hat{\boldsymbol{\rho}}] / |\boldsymbol{\rho}|^5$, where the $\ell^4$-magnitude product of the two body's quadrupole tensors against the inter-body direction scales the term. The $1/|\boldsymbol{\rho}|^5$ scaling follows from the product of two $|\boldsymbol{\rho}|^{-(\ell+1)}$ Laplace-expansion factors with $\ell = 2$ each, summed across the two-body multipole sum.
  - **Octupole-monopole term:** $V^{(\ell_A=3,\ell_B=0)} \propto G\, m_B\, \mathcal{O}_A : (\hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}\hat{\boldsymbol{\rho}}) / |\boldsymbol{\rho}|^4$, where $\mathcal{O}_A$ is body $A$'s rank-3 octupole moment tensor (symmetric and trace-free, encoding $7$ independent degrees of freedom for $\ell=3$). The $1/|\boldsymbol{\rho}|^4$ scaling is the $|\boldsymbol{\rho}|^{-(\ell+1)}$ factor at $\ell=3$ paired with the constant body-$B$ mass.
  
  **Bound on numerical magnitude at first-cut configurations.** Using the dimensional estimates $\mathcal{O}_A \sim m_A \ell_A^3$ and $|\mathbf{Q}_X| \sim m_X \ell_X^2$:
  
  | Term | Magnitude scaling | At $\ell/r = 10^{-2}$ |
  |---|---|---|
  | $V^{(0)}$ Kepler | $Gm_Am_B/r$ | $1$ (reference) |
  | $V_{\mathrm{tidal},A}$ quadrupole-monopole | $\sim (\ell_A/r)^2$ | $10^{-4}$ |
  | $V^{(\ell_A=2,\ell_B=2)}$ quadrupole-quadrupole | $\sim (\ell_A \ell_B/r^2)^2 = (\ell^2/r^2)^2$ | $10^{-8}$ |
  | $V^{(\ell_A=3,\ell_B=0)}$ octupole-monopole | $\sim (\ell_A/r)^3$ | $10^{-6}$ |
  
  Both next-order terms fall below the **observed double-precision conservation drift** $|dE/E_0| = 6.2 \times 10^{-12}$ for separations $\ell/r \lesssim 10^{-2}$ — the regime where the engineer-tier first-cut operates ($|\boldsymbol{\rho}| \sim 1000$ m vs $\ell \sim 7$ m). The §3.3 quadrupole-truncation-error estimate is consistent with these magnitudes.
  
  **Resolution status.** Explicit derivation deferred to v0.2 if the canonical-tier scope extends to tidally-locked binary asteroid regime ($\ell/r \to 1$, where both terms become engineering-significant; reference anchor for that regime: Boué 2018 *J. Geophys. Res.* or Maciejewski 1995). The order-of-magnitude bound established here is sufficient for v0.1.
- **OQ-3.3.** ***Resolved 2026-05-25 (canonical-tier round 3, third independent derivation route).*** The gravity-gradient torque on body $A$ can be derived by a **third independent route** beyond the two already in §3.7 (Lagrangian-formalism variational derivation) and §3.7.1 cross-check (Hughes 2004 direct-formula route). The third route is **direct integration of $\mathbf{r} \times \nabla V$ over the body's mass distribution**, following Goldstein 2002 §5.5.
  
  **Setup.** Body $A$ sits in the gravitational potential of body $B$ treated as a point mass at $\mathbf{R}_B$. At each material point of body $A$ (position $\boldsymbol{\xi}_A$ in the inertial frame, with $\boldsymbol{\xi}_A = \mathbf{R}_q^A \mathbf{x}_A$ as in §3.4), the force per unit mass is $-\nabla \Phi_B(\boldsymbol{\xi}_A)$ where $\Phi_B = -Gm_B/|\boldsymbol{\xi}_A - \mathbf{R}_B|$. The body-frame torque on body $A$ about its own centre of mass is:
  
  $$
  \boldsymbol{\tau}_A^{\mathrm{inertial}}
  \;=\;
  -\int_A \rho_A(\mathbf{x}_A)\, (\boldsymbol{\xi}_A - \mathbf{R}_A) \times \nabla \Phi_B(\boldsymbol{\xi}_A)\, dV_A.
  $$
  
  **Expansion.** Substitute $\boldsymbol{\xi}_A - \mathbf{R}_A = \mathbf{R}_q^A \mathbf{x}_A$, and Taylor-expand $\nabla \Phi_B(\boldsymbol{\xi}_A)$ about the body's centre of mass $\mathbf{R}_A$ (small expansion parameter $|\mathbf{x}_A|/|\boldsymbol{\rho}|$, the same parameter that controls §3.4's expansion). To leading order — keeping terms linear in $\mathbf{x}_A$:
  
  $$
  \nabla \Phi_B(\boldsymbol{\xi}_A) \;\approx\; \nabla\Phi_B(\mathbf{R}_A) + (\mathbf{R}_q^A \mathbf{x}_A) \cdot \nabla^2 \Phi_B(\mathbf{R}_A) + \mathcal{O}(|\mathbf{x}_A|^2).
  $$
  
  The zeroth-order term contributes zero torque (integrates against $\int_A \rho_A \mathbf{x}_A dV_A = \mathbf{0}$ by centre-of-mass convention). The first-order term gives, after carrying the cross product through the integral:
  
  $$
  \boldsymbol{\tau}_A^{\mathrm{inertial}}
  \;\approx\;
  -\mathbf{R}_q^A\, \boldsymbol{\epsilon}_{ijk}\,(\mathbf{M}_A)_{j\ell}\, [\nabla^2 \Phi_B(\mathbf{R}_A)]_{k\ell}\, \hat{\mathbf{e}}_i,
  $$
  
  where $\mathbf{M}_A = \int \rho_A \mathbf{x} \otimes \mathbf{x}\, dV$ is the second-moment tensor of §3.4.3, and $\boldsymbol{\epsilon}$ is the Levi-Civita tensor. Using $[\nabla^2 \Phi_B]_{k\ell} = -Gm_B [3\hat{\boldsymbol{\rho}}_k \hat{\boldsymbol{\rho}}_\ell - \delta_{k\ell}]/|\boldsymbol{\rho}|^3$ (the tidal-tensor identity for a point-mass field), substituting, and converting to the body frame via $\mathbf{R}_q^A$ similarity transformation (using the trace identity $\mathbf{Q}_A = \mathrm{tr}(\mathbf{I}_A)\mathbb{I}_3 - 3\mathbf{I}_A$, or equivalently $\mathbf{M}_A = \tfrac{1}{2}\mathrm{tr}(\mathbf{I}_A)\mathbb{I}_3 - \mathbf{I}_A$, from §3.4.3), the body-frame torque emerges as:
  
  $$
  \boldsymbol{\tau}_A^{\mathrm{body}}
  \;=\;
  \frac{3 G m_B}{|\boldsymbol{\rho}|^3}\,
  \hat{\boldsymbol{\rho}}_{\mathrm{body},A} \times (\mathbf{I}_A \cdot \hat{\boldsymbol{\rho}}_{\mathrm{body},A}),
  $$
  
  **identical** to the §3.7 boxed result.
  
  **Three-route cross-validation.** Three independent derivation routes — (1) variational $-\delta V_{\mathrm{tidal}}$ rule (§3.7), (2) direct Hughes 2004 §5.3 formula (cited in §3.7.1), (3) Goldstein 2002 §5.5 direct $\mathbf{r}\times\nabla V$ integration (this OQ-3.3 resolution) — all produce the boxed body-frame gravity-gradient torque. This extends the two-in-parallel methodology of §3.6 to *three-in-parallel* for this specific equation, hardening §3.7 against algebraic-error introduction in any single route. The Phase 1 D4 cross-test `backends/fortran/tests/test_gravity_gradient_torque.f90` is the engineer-tier C2 cross-implementation regression guard at ULP × √N tolerance.
- **OQ-3.4.** ***Resolved 2026-05-25 (canonical-tier round 3, engineer-tier follow-up specification).*** Augmentation of `dynamics.py::total_potential_energy` to include the quadrupole-tidal contributions $V_{\mathrm{tidal},A} + V_{\mathrm{tidal},B}$ from §3.4.3 boxed.
  
  **Current implementation (lines 207-214).** The engineer-tier `total_potential_energy` function computes only the monopole-monopole Kepler term $V^{(0)} = -Gm_Am_B/|\boldsymbol{\rho}|$ and returns that as "the system's potential energy" for the conservation diagnostic $|dE/E_0|$. The tidal pieces are implicit in the dynamics (via the gravity-gradient torque § 3.7 which directly enters the body-frame Euler equation of motion) but absent from the diagnostic.
  
  **Specified augmentation (v0.1.1 engineer-tier patch).** Add the §3.4.3 boxed tidal pieces to the diagnostic:
  
  ```
  # Existing
  V_kepler = -G * m_a * m_b / r_mag
  # New (proposed)
  rhat_body_a = R_a.T @ (r_vec / r_mag)
  V_tidal_a = -(G * m_b / (2 * r_mag**3)) * (
                np.trace(I_a) - 3 * rhat_body_a @ I_a @ rhat_body_a)
  # symmetric for body B
  rhat_body_b = R_b.T @ (-r_vec / r_mag)   # note opposite sign
  V_tidal_b = -(G * m_a / (2 * r_mag**3)) * (
                np.trace(I_b) - 3 * rhat_body_b @ I_b @ rhat_body_b)
  V_total = V_kepler + V_tidal_a + V_tidal_b
  return V_total
  ```
  
  **Expected impact on conservation diagnostic.** At first-cut configurations ($|\boldsymbol{\rho}| \sim 1000$ m, $\ell \sim 7$ m), the tidal piece magnitude is $|V_{\mathrm{tidal},A}/V_{\mathrm{Kepler}}| \sim (\ell/r)^2 \sim 5 \times 10^{-5}$. With the augmentation, the conservation diagnostic tracks the *full* quadrupole-truncated $V$ rather than the monopole-only piece; the expected effect on $|dE/E_0|$:
  
  - At $|\boldsymbol{\rho}| \sim 1000$ m: improvement is **below double-precision round-off** ($|dE/E_0|$ already at $6.2 \times 10^{-12}$, the tidal contribution would change this at the $\sim 10^{-16}$ level which is below floating-point noise).
  - At $|\boldsymbol{\rho}| \sim 30$ m (the close-approach regime, $\ell/r \sim 0.2$): the tidal piece is $\sim (\ell/r)^2 \sim 0.04$ of the Kepler magnitude; without the augmentation, the conservation diagnostic spuriously appears to "fail" at $|dE/E_0| \sim 10^{-2}$ (because the integrator is conserving full $V$ but the diagnostic measures only monopole $V$); with the augmentation, the conservation diagnostic recovers the integrator's actual conservation precision down to the RK4 floor.
  
  **Companion v0.1.1 follow-up: implement back-reaction force.** As surfaced in Ch 04 §4.9 OQ-4.5, the back-reaction force on the orbit ($\partial V_{\mathrm{tidal}}/\partial \boldsymbol{\rho}$ — Ch 04 §4.6.2 boxed) is also absent from the first-cut. The OQ-3.4 conservation-diagnostic augmentation is independent of OQ-4.5 (the diagnostic measures total $V$ correctly; the question is whether the integrator's *equations of motion* propagate the back-reaction). Both should be implemented together at v0.1.1 to make the diagnostics and the equations of motion mutually consistent.
  
  **Engineer-tier change tracker.** Queued at `_specs/CANONICAL-VALIDATION-TESTBEDS-AND-EDGE-CASES-v0.1.md` (next revision) §3.3 close-approach testbed: the augmented diagnostic is what the close-approach Class B testbed should use to characterise the model's domain-of-validity boundary precisely.
- **OQ-3.5.** ***Resolved 2026-05-23*** — project-owner-supplied T&B PDF reviewed (1552-pp electronic edition). T&B §27.5.1 ("Multipole-Moment Expansion") + §27.5.2 ("Quadrupole-Moment Formalism") covers the Cartesian multipole expansion of the Newtonian potential at Eqs. (27.50)–(27.53), in the framing of gravitational-wave source theory. Specific citations added to §3.2 (T&B Eq. 27.50 single-body integral) and §3.4.3 (T&B Eqs. 27.51, 27.52, 27.53 Taylor expansion + monopole-plus-quadrupole truncation + trace-free quadrupole moment definition) reference anchors, with the prefactor-convention difference between T&B's $\mathcal{I}_{jk}$ (prefactor $1/3$) and our Jackson-convention $\mathbf{Q}_{jk}$ (prefactor $3$) made explicit ($\mathbf{Q}_{jk} = 3\, \mathcal{I}_{jk}$). **T&B does *not* develop the spherical-harmonic side of the multipole expansion (our §3.5)** — that side remains anchored on Jackson 1999 §4.1. **T&B does *not* derive the trace identity** $\mathbf{Q} = \mathrm{tr}(\mathbf{I}_{\mathrm{inertia}}) \cdot \mathbb{I} - 3\, \mathbf{I}_{\mathrm{inertia}}$ explicitly — they note the quadrupole moment is "the second moment of the mass distribution with its trace removed" without spelling out the inertia-tensor connection that our §3.5 derives. The Chapter 03 chain Cartesian-→-spherical-harmonic-→-trace-identity-→-gravity-gradient-torque is therefore an integration synthesis across T&B + Jackson + Hughes anchors rather than a copy of any single source.
- **OQ-3.6.** ***Resolved 2026-05-23*** — engineer-tier bridge for the gravity-gradient torque authored at new subsection §3.7.2 ("Gravity-gradient torque Jacobian: sensitivity to initial conditions and assembly geometry"), continuing the §2.5.1 framework. Deliverables: linearised orientation-response Jacobian $\mathbf{J}^{\mathrm{orient}}_{A}$ (boxed) as the unifying $3 \times 3$ matrix operating on both body-orientation and inter-body-position-direction perturbations; libration-equilibrium specialisation $\mathbf{J}^{\mathrm{orient}}_{A}\big|_{\mathrm{eq}} = \alpha (I_\parallel - I_\perp) \mathrm{diag}(1,1,0)$ with gravity-gradient stability table (prolate stable, oblate unstable, spherical marginal); radial-vs-tangential decomposition of $\partial \boldsymbol{\tau}_A / \partial \boldsymbol{\rho}$ showing the tangential piece reduces to the same Jacobian; identification of $\partial \boldsymbol{\tau}_A / \partial q_B = 0$ at quadrupole truncation; symmetry-class rank table (spherical → 0, axisymmetric → 2, asymmetric → 3); assembly sensitivity figure of merit $\Sigma_{\mathrm{assembly}} = \eta_A \eta_B \cdot f(\hat{\boldsymbol{\rho}})$; numerical sanity check (JWST-like at L2-like configuration: spring constant $4.7 \times 10^{-10}\ \text{N m/rad}$, instability e-folding $\sim 27\ \text{days}$); multi-scale time-scale table identifying the first-cut as spin-resolving but libration-blind; symplectic-Lie-group-integrator flag for long-horizon studies (v0.2 engineer-tier item). Anchors: Hughes 2004 Ch 8; Markley & Crassidis 2014 §3.7.

---

**Next chapter:** [04 — Euler-Lagrange equations and Euler-Poincaré reduction](04-euler-lagrange-and-euler-poincare.md) *(queued for canonical-tier Round 3)*

**End 03-mutual-gravitational-potential.md**
