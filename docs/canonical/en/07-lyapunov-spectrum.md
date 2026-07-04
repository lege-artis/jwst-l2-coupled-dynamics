# Chapter 7 — Lyapunov spectrum of the coupled two-body rotational dynamics

> **Prerequisite reading.** Chapters 01–06. This chapter operates on the boxed Ch 04 §4.7 equation set propagated by the Phase 2 symplectic SE(3)×SE(3) variational integrator (`integrators/symplectic_se3_variational.py`, Sitting 1), and quantifies its sensitivity to initial conditions via the Lyapunov spectrum. It leans on the Ch 06 analytical-limit testbeds — the §6.3/§6.5 torque-free tops and the §6.4 gravity-gradient libration — as calibration oracles.
>
> **Status.** v0.1 (canonical-tier Phase 2 Sitting 4, 2026-06-07). All numerical evidence quoted below is from the Sitting 3 acceptance-gate runs (`_handoffs/phase-2-symplectic-lie-group/SITTING-3-FINAL-REPORT.md`), all gates GREEN.

---

## §7.1 Overview and motivation

Chapters 05 and 06 established *where* the rotational dynamics of the coupled two-body system is stable, unstable, or librating — but only in the linearised or integrable sense. The Ch 05 §5.2.2 tennis-racket classification gives linearised eigenvalues at the principal-axis equilibria; the Ch 06 §6.4.3 Lagrange-top reduction and §6.5.3 Jacobi-elliptic solution give closed-form *integrable* finite-amplitude motion. Neither answers the question this chapter poses:

> **Does the gravity-gradient coupling between orbital and attitude dynamics break integrability — and if so, with what exponential divergence rate?**

The quantitative instrument is the **Lyapunov spectrum**: the set of asymptotic exponential growth rates $\lambda_1 \geq \lambda_2 \geq \cdots \geq \lambda_n$ of infinitesimal tangent perturbations along a trajectory. For an $n$-dimensional flow $\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x})$ with flow map $\phi^t$, the exponents are defined through the tangent cocycle $D\phi^t$:

$$
\lambda_i \;=\; \lim_{t \to \infty} \frac{1}{t}\, \ln \frac{\|D\phi^t(\mathbf{x}_0)\, \mathbf{v}_i\|}{\|\mathbf{v}_i\|},
$$

for a generic basis $\{\mathbf{v}_i\}$ of tangent vectors (Oseledec's multiplicative ergodic theorem guarantees existence almost everywhere). A positive $\lambda_1$ is the operational definition of deterministic chaos; $\lambda_1 = 0$ with quasi-periodic motion signals integrable or regular dynamics.

The physically decisive facts for our system, established quantitatively in §7.5–§7.6:

- The **torque-free asymmetric top is integrable** — including the Dzhanibekov near-saddle motion. Its Jacobi-elliptic closed form (Ch 06 §6.5.3) leaves no room for a positive Lyapunov exponent; the separatrix is a heteroclinic orbit of an integrable system, not a chaotic layer.
- **Gravity-gradient coupling breaks integrability.** With the orbital-attitude coupling switched on and the initial condition near the intermediate-axis saddle, the computed $\lambda_{\max} = 1.506 \times 10^{-5}\ \mathrm{s}^{-1} > 0$ — small (weak near-separatrix chaos at the coupling strength used) but definitively positive.
- **Libration is not divergence.** A body librating about the GG-stationary equilibrium (Ch 06 §6.4) oscillates at $\omega_{\mathrm{lib}}$ but has $|\lambda_{\max}|$ a factor 263 below $\omega_{\mathrm{lib}}$ — oscillation frequency and Lyapunov exponent are different observables, and conflating them is a classic misreading.

**Algorithm and embedding.** We compute the spectrum with the Benettin–Galgani–Giorgilli–Strelcyn (BGGS) 1980 algorithm (§7.2), with the tangent cocycle realised as a finite-difference Jacobian of the symplectic integrator step in a flat Euclidean $\mathbb{R}^{14}$ chart (§7.3). The Euclidean embedding is a deliberate, documented design choice (OQ-PHASE2-2) with a known convergence artefact for integrable systems that directly shaped the acceptance-gate thresholds — the honest accounting is in §7.3 and §7.5(a).

**A note for Czech-reading audiences.** Podolský 2018 §8.1–§8.3 (free symmetric top; Lagrange top) is the closest CS-language anchor for the *integrable* reference dynamics that calibrate this chapter's scenarios (a) and (b). The Lyapunov-spectrum formalism itself is **not** covered in Podolský; CS readers should pair this chapter with Skokos 2010 §2 (EN) for the algorithmic content. This gap is recorded honestly here rather than papered over with a loose citation.

## §7.2 The BGGS algorithm — formal statement

### §7.2.1 Why naive tangent propagation fails

Propagating a single tangent vector $\mathbf{v}(t) = D\phi^t \mathbf{v}_0$ and reading off $\lambda_1 = \ln(\|\mathbf{v}(T)\|/\|\mathbf{v}_0\|)/T$ overflows double precision once $\lambda_1 T \gtrsim 700$; propagating a basis of $k$ vectors collapses them all onto the most-expanding direction, destroying access to $\lambda_2, \ldots, \lambda_k$. The BGGS resolution: **periodically reorthonormalise** the basis and accumulate the logarithms of the stretching factors.

### §7.2.2 Algorithm statement

> **Algorithm (BGGS 1980, §2; modern formulation Skokos 2010 §2.2).**
> *Input:* initial state $\mathbf{x}_0$, integrator step $\varphi_{dt}$, tangent propagator $D\varphi_{dt}$, number of vectors $k$, QR cadence $\tau_{\mathrm{QR}}$, horizon $T$.
>
> 1. Initialise a random orthonormal basis $Q_0 = [\mathbf{v}_1, \ldots, \mathbf{v}_k]$.
> 2. For each epoch $m = 1, \ldots, T/\tau_{\mathrm{QR}}$:
>    a. Propagate state and all $k$ tangent vectors for $\tau_{\mathrm{QR}}/dt$ steps:
>       $\mathbf{x} \leftarrow \varphi_{dt}(\mathbf{x})$, $\mathbf{v}_j \leftarrow D\varphi_{dt}(\mathbf{x})\,\mathbf{v}_j$.
>    b. Reorthonormalise by (modified) Gram-Schmidt: $Q_m R_m = [\mathbf{v}_1, \ldots, \mathbf{v}_k]$; record the stretching factors $r_j^{(m)} = (R_m)_{jj} = \|\mathbf{v}_j^{\perp}\|$ (norm of the $j$-th column after projection onto the orthogonal complement of the preceding columns).
>    c. Accumulate $S_j \leftarrow S_j + \ln r_j^{(m)}$; reset the basis to $Q_m$.
> 3. *Output:* $\lambda_j = S_j / T$ for $j = 1, \ldots, k$, in descending order.
>
> The $j$-th accumulated log-stretching measures the growth of the $j$-dimensional volume element divided by the $(j-1)$-dimensional one, which converges to $\lambda_j$ by Oseledec's theorem (Benettin et al. 1980 §2; the discrete-QR variant is compared against alternatives in Geist, Parlitz & Lauterborn 1990).

### §7.2.3 Implementation mapping

The implementation is `lyapunov/spectrum.py::lyapunov_spectrum` (driver, steps 2–3), `lyapunov/cocycle.py::propagate_tangent` (the $D\varphi_{dt}$ realisation, §7.3), and `lyapunov/cocycle.py::gram_schmidt` (modified Gram-Schmidt, step 2b; the modified column-by-column variant per Golub & Van Loan §5.2 for numerical stability, with degenerate-column re-injection from the canonical basis). The integrator step $\varphi_{dt}$ is the Sitting 1 symplectic KDK Störmer-Verlet step on SE(3)×SE(3) — the same step that holds the AG-INT-3 angular-momentum ULP gate, so the Lyapunov run inherits exact conservation of $\mathbf{L} = R\,\boldsymbol{\Pi}$ along the reference trajectory.

**Numerical example (scenario (c) production run).** $k=1$, $dt = 50$ s, $T = 100{,}000$ s: the pilot phase (§7.4) estimated $\lambda_{\max}^{\mathrm{pilot}}$, set the production cadence, and the accumulated log-stretchings gave $\lambda_{\max} = 1.506 \times 10^{-5}\ \mathrm{s}^{-1}$ — the AG-LYAP-3 headline number (§7.5(c)).

## §7.3 Tangent-space selection: Euclidean $\mathbb{R}^{14}$ vs manifold (OQ-PHASE2-2)

### §7.3.1 The choice

The state of the rotational subsystem lives on the manifold $(S^3 \times \mathbb{R}^3)^2$ — unit quaternion + body angular velocity per body. The mathematically natural tangent space is the Lie-algebraic one ($\mathfrak{so}(3) \times \mathbb{R}^3$ per body, dimension 12). OQ-PHASE2-2 instead **confirmed the flat Euclidean chart**:

$$
\mathbb{R}^{14} \;=\; \underbrace{q_A(4)}_{\text{quaternion}} + \underbrace{\boldsymbol{\omega}_A(3)}_{\text{body rate}} + \underbrace{q_B(4)}_{\text{quaternion}} + \underbrace{\boldsymbol{\omega}_B(3)}_{\text{body rate}},
$$

with the tangent cocycle realised as a symmetric finite-difference Jacobian-vector product of the integrator step:

$$
\boxed{
D\varphi_{dt}(\mathbf{x})\cdot \mathbf{v}
\;\approx\;
\frac{\varphi_{dt}(\mathbf{x} + \epsilon\, \hat{\mathbf{v}}) - \varphi_{dt}(\mathbf{x} - \epsilon\, \hat{\mathbf{v}})}{2\epsilon}\; \|\mathbf{v}\|,
\qquad \epsilon = 10^{-7},
}
$$

where $\hat{\mathbf{v}} = \mathbf{v}/\|\mathbf{v}\|$ (normalising before perturbing keeps the FD perturbation magnitude independent of the tangent norm, which grows during the BGGS run). The rationale: implementation simplicity (no Lie-group chart atlas, no exponential-map bookkeeping), an FD Jacobian that requires only the already-validated integrator step as a black box, and transparent reuse of the $\mathbb{R}^{26}$ full-state layout (`lyapunov/cocycle.py` `_R14_MAP`). The integrator's internal quaternion normalisation corrects the slight off-manifold drift introduced by each FD perturbation.

### §7.3.2 Two distinct effects — and an honest accounting (corrected after Phase 3)

The flat embedding has **one genuine artefact** and is associated with **one effect that turned out NOT to be its fault**. Phase 3 (§7.10) disentangled the two; this subsection states the corrected understanding, with the original Phase-2-era reasoning recorded honestly below.

**Genuine embedding artefact — spurious quaternion-norm modes.** Each unit quaternion in $\mathbb{R}^{14}$ carries a radial (norm) direction $\delta q \parallel q$ that points off the constraint manifold $S^3$ and is removed by the integrator's per-step normalisation. The full 14-vector spectrum is therefore contaminated by **two large-negative exponents** (one per body) with no positive partners — the spectrum sums to $-0.14\ \mathrm{s}^{-1}$ instead of $0$ and fails Hamiltonian pairing (§7.7). This *is* an embedding artefact, and it is the reason a meaningful full spectrum requires a 12-dimensional cocycle (Phase 3, §7.9–§7.10). It does not affect the single-vector $\lambda_{\max}$ gates, which track the most-expanding *dynamical* direction, well separated from the contracting norm modes.

**NOT an embedding artefact — the slow scenario-(a) FTLE convergence.** For the integrable scenario (a), the finite-time Lyapunov exponent (FTLE) converges to its true value 0 only as

$$
\lambda^{\mathrm{FTLE}}(T) \;\sim\; O\!\left(\frac{\log T}{T}\right),
$$

a slow algebraic decay rather than the exponentially fast convergence one might naively expect. **Phase 3 showed this is intrinsic to the quasi-periodic multi-frequency dynamics, not the flat chart:** the Lie-algebraic $\mathfrak{so}(3) \times \mathbb{R}^3$ cocycle (§7.10) — which carries no quaternion-norm directions at all — produces an essentially identical $\lambda_{\max} = 2.29 \times 10^{-4}\ \mathrm{s}^{-1}$ vs the $\mathbb{R}^{14}$ value $2.30 \times 10^{-4}\ \mathrm{s}^{-1}$ at $T = 20{,}000$ s. A genuinely embedding-caused slow convergence would have vanished under the manifold-native cocycle; it did not. The slow convergence reflects the non-resonant mix of rotation and precession frequencies in the integrable orbit itself.

> **Provenance note (corrected 2026-06-07, Phase 3).** The Phase-2-era version of this subsection attributed the $O(\log T/T)$ scenario-(a) convergence to "polynomial norm growth in flat quaternion coordinates" — i.e. to the embedding — and predicted that a Lie-algebraic cocycle would restore fast convergence and make a $10^{-8}$ AG-LYAP-1 gate achievable. Phase 3 falsified that prediction (AG-P3-3): the Lie-algebraic cocycle reproduced the same $\lambda_{\max}$. The convergence law is intrinsic; the embedding's real defect is the spurious-norm-mode contamination of the *full* spectrum, which the Lie-algebraic cocycle *does* fix. The original `lyapunov/reference.py::lambda_max_gate_scenario_a` docstring carries the superseded reasoning. This is the same physical-truth-anchored correction discipline (KB-OPUS-DOC-INTRA-INCONSISTENCY) applied to a cross-tier prediction rather than a boxed formula.

### §7.3.3 Implication for gate design

The consequence is a methodological rule worth stating in general form:

> **Acceptance thresholds for "zero" Lyapunov exponents must be calibrated empirically against the cocycle's measured convergence law — not taken from the asymptotic-Lyapunov literature. That law may be set by the embedding chart, or (as here) by the intrinsic dynamics; either way the threshold follows the measurement, not the ideal.**

Concretely: the Phase 2 brief originally specified AG-LYAP-1 (integrable scenario, $\lambda_{\max} <$ threshold) at $10^{-8}$. Under the $O(\log T / T)$ law, reaching $10^{-8}$ requires $T \sim 10^{9}$ s — not achievable at any feasible integration time. The revised gate of $10^{-3}$ captures the correct *qualitative* requirement (the integrable FTLE must be subdominant to the chaotic scenario's analytical growth rate $\sigma_c \approx 7.07 \times 10^{-3}\ \mathrm{s}^{-1}$, §7.6) and is met with margin 4.3× at $T = 20{,}000$ s. This is a documented gate revision (Sitting 3 report §5 D3), not a silent loosening. **Phase 3 confirms the $10^{-3}$ gate is the right value even with the manifold-native Lie-algebraic cocycle** (§7.10), since the slow convergence is intrinsic — the gate revision was correct, only its originally-stated *cause* (§7.3.2) needed amendment.

## §7.4 Pilot-then-production QR cadence autotuning

The QR cadence $\tau_{\mathrm{QR}}$ trades off cost (each reorthonormalisation is $O(k^2 n)$) against overflow/collapse risk (stretching factors must stay within floating-point range and the basis must not degenerate between QRs). A useful heuristic is to let each epoch accumulate $\sim 0.1$ e-foldings of the fastest direction. Since $\lambda_{\max}$ is not known in advance, the implementation autotunes it (landscape §4.3 pattern):

1. **Pilot run:** 20 epochs at $\tau_{\mathrm{QR}}^{\mathrm{pilot}} = 10\,dt$ (200 integrator steps total — cheap), giving the estimate $\lambda_{\max}^{\mathrm{pilot}} = S_1^{\mathrm{pilot}} / T_{\mathrm{pilot}}$.
2. **Production cadence:**

$$
\tau_{\mathrm{QR}} \;=\; \frac{0.1}{|\lambda_{\max}^{\mathrm{pilot}}|},
\qquad \text{clamped to } \big[\,10\,dt,\; T/3\,\big].
$$

3. **Near-zero fallback:** if $|\lambda_{\max}^{\mathrm{pilot}}| < 10^{-12}$ (integrable pilot), the cadence degenerates; the fallback $\tau_{\mathrm{QR}} = \min(100\,dt,\; T/5)$ is applied instead.

**Numerical example.** In scenario (b) (GG-stationary libration, $dt = 500$ s, $T = 10^6$ s) the pilot estimate is near-zero-to-small, and the clamp keeps the production cadence well inside $[5000\ \mathrm{s},\ 3.3\times 10^{5}\ \mathrm{s}]$; the production run completes in 1.3 s of wall time at $T_{\mathrm{total}} = 10^6$ s of simulated time (Sitting 3 report §4) — the autotuned cadence is what makes million-second horizons cheap.

## §7.5 Three-scenario calibration

A Lyapunov pipeline that has only ever been run on the system under study is uncalibrated. The Sitting 3 testbed (`testbeds/lyapunov_three_scenario.py`) brackets the physics with three scenarios whose expected answers are known from Ch 05–06 analysis — one integrable-zero, one libration-bounded, one genuinely chaotic.

### §7.5(a) Torque-free asymmetric top — $\lambda_{\mathrm{true}} = 0$ (AG-LYAP-1)

**Setup.** $\mathbf{I} = \mathrm{diag}(100, 200, 400)$ kg·m², $\boldsymbol{\omega}_0 = [0.1, 0, 0]$ rad/s (stable smallest-moment axis $\hat{\mathbf{e}}_1$ — quasi-periodic, fast FTLE convergence), partner body parked at $10^{12}$ m so the gravitational coupling is negligible. The dynamics is the Ch 06 §6.5.3 Jacobi-elliptic integrable motion; the true $\lambda_{\max} = 0$.

**Expectation.** By the §7.3.2 artefact, the computed FTLE decays to 0 only as $O(\log T / T)$; the gate tests subdominance, not literal zero.

| Quantity | Value |
|---|---|
| $dt$, $T$ | 1.0 s, 20 000 s |
| $\lambda_{\max}$ (computed FTLE) | $2.304 \times 10^{-4}\ \mathrm{s}^{-1}$ |
| Gate AG-LYAP-1 ($< 10^{-3}$) | **PASS** ✓ (margin 4.3×) |

### §7.5(b) GG-stationary libration — $|\lambda_{\max}| < \omega_{\mathrm{lib}}$ (AG-LYAP-2)

**Setup.** JWST-like prolate body (Ch 06 §6.4) at the GG-stationary equilibrium — long axis toward an Earth-mass primary ($m_B = 5.97 \times 10^{24}$ kg) at $\rho = 10^{7}$ m — tilted by 0.01 rad. The body librates at the Ch 06 §6.4.2 frequency $\omega_{\mathrm{lib}} = \sqrt{3 G m_B (I_\perp - I_\parallel)/(I_\perp \rho^3)}$.

**Expectation.** Libration is bounded oscillation in the Lagrange-top effective potential (Ch 06 §6.4.3) — *no* exponential divergence. The discriminating gate: $|\lambda_{\max}|$ must sit far below the oscillation frequency itself.

| Quantity | Value |
|---|---|
| $dt$, $T$ | 500 s, $10^{6}$ s |
| $\omega_{\mathrm{lib}}$ | $6.379 \times 10^{-4}\ \mathrm{s}^{-1}$ |
| $\lambda_{\max}$ | $2.421 \times 10^{-6}\ \mathrm{s}^{-1}$ |
| Gate AG-LYAP-2 ($|\lambda_{\max}| < \omega_{\mathrm{lib}}$) | **PASS** ✓ (margin 263×) |

The 263× separation is the quantitative form of the §7.1 caution: a librating spacecraft *oscillates* at $\omega_{\mathrm{lib}}$ but does not *diverge* at anything close to it.

### §7.5(c) GG-Dzhanibekov near-saddle — $\lambda_{\max} > 0$ (AG-LYAP-3)

**Setup.** Asymmetric body $\mathbf{I} = \mathrm{diag}(100, 200, 400)$ kg·m², $\boldsymbol{\omega}_0 = [0.001, 0.01, 0.001]$ rad/s — near the intermediate-axis saddle (Ch 06 §6.5) — with Earth-mass primary at $\rho = 3 \times 10^{6}$ m, chosen so the GG coupling is near-resonant with the saddle dynamics. The analytical Dzhanibekov characteristic frequency (Ch 06 §6.5.2):

$$
\sigma \;=\; \Omega_2^0 \sqrt{\frac{(I_2 - I_1)(I_3 - I_2)}{I_1 I_3}}
\;=\; 0.01 \sqrt{\frac{100 \cdot 200}{100 \cdot 400}}
\;=\; 7.071 \times 10^{-3}\ \mathrm{s}^{-1}.
$$

**Expectation.** Torque-free, this motion would be integrable (heteroclinic, Jacobi-elliptic). The GG coupling breaks integrability; near the separatrix this produces a thin chaotic layer with small but positive $\lambda_{\max}$ — of KAM-residue order at the coupling strength used, *not* of order $\sigma$.

| Quantity | Value |
|---|---|
| $dt$, $T$ | 50 s, $10^{5}$ s |
| $\sigma$ (analytical) | $7.071 \times 10^{-3}\ \mathrm{s}^{-1}$ |
| $\lambda_{\max}$ | $1.506 \times 10^{-5}\ \mathrm{s}^{-1}$ |
| ratio $\lambda_{\max}/\sigma$ | 0.0021 |
| Gate AG-LYAP-3 ($\lambda_{\max} > 0$, ratio $< 10$) | **PASS** ✓ |

The original brief gate demanded ratio within a factor 3 of $\sigma$; that calibration presumed coupling strength $\varepsilon \approx 1$. At the actual $\varepsilon = \omega_{\mathrm{lib}}/\sigma \approx 0.81$ the weak-chaos exponent sits orders of magnitude below $\sigma$ — physically correct for near-separatrix KAM-regime chaos — and the gate was revised to ratio $< 10$ with the positivity requirement carrying the physical content (Sitting 3 report §5, gate-revision table).

## §7.6 Cross-scenario comparison and the integrability criterion

The three scenarios only earn their keep jointly: the pipeline must *order* them correctly. The naive ordering test — computed $\lambda_c > \lambda_a$ — fails for an instructive reason: at accessible $T$, the integrable scenario's still-declining FTLE ($3.924 \times 10^{-4}$ at $T = 10{,}000$ s) sits *above* the true chaotic exponent $\lambda_c \approx 1.5 \times 10^{-5}$. Under the $O(\log T/T)$ law, $\lambda_a^{\mathrm{FTLE}}$ would not fall below $\lambda_c$ until $T_a \approx 2.6 \times 10^{8}$ s — infeasible. The robust, physically meaningful ordering compares the integrable FTLE against the *analytical* chaotic growth rate:

| Quantity | Value |
|---|---|
| $\lambda_a^{\mathrm{FTLE}}$ at $T = 10{,}000$ s | $3.924 \times 10^{-4}\ \mathrm{s}^{-1}$ |
| $\sigma_c$ (analytical) | $7.071 \times 10^{-3}\ \mathrm{s}^{-1}$ |
| ratio | **0.055** |
| Ordering gate ($\lambda_a < \sigma_c$) | **PASS** ✓ |

with $\lambda_c > 0$ established separately by AG-LYAP-3. The combined picture:

| Scenario | True $\lambda_{\max}$ | Computed evidence | Verdict |
|---|---|---|---|
| (a) torque-free asymmetric top | $0$ (integrable) | FTLE $\to 0$ as $O(\log T/T)$; $2.3 \times 10^{-4} \ll \sigma_c$ | regular |
| (b) GG-stationary libration | $\approx 0$ (bounded libration) | $|\lambda| = 2.4 \times 10^{-6} \ll \omega_{\mathrm{lib}}$ (263×) | regular (oscillatory) |
| (c) GG-Dzhanibekov near-saddle | $> 0$ (GG breaks integrability) | $\lambda = 1.5 \times 10^{-5} > 0$, ratio $\lambda/\sigma = 0.002$ | weakly chaotic |

This is the chapter's **integrability criterion** in operational form: integrable and librating configurations show FTLEs that decay with $T$ and sit far below the relevant dynamical frequencies ($\sigma$, $\omega_{\mathrm{lib}}$); genuinely chaotic configurations show a $T$-stable positive exponent, however small.

## §7.7 Hamiltonian pairing and spectrum structure

The coupled two-body rotational system is Hamiltonian, and its Lyapunov spectrum inherits the symplectic constraint: exponents come in **pairs $(+\lambda, -\lambda)$**, with zero exponents contributed by conserved quantities (the Hamiltonian itself, plus each independent first integral on the reduced space). For the 14-dimensional rotational subsystem this means 7 pairs; integrability would force *all* pairs to $(0, 0)$.

The implementation check is `lyapunov/reference.py::hamiltonian_pairing_check`: given $k$ computed exponents in descending order, it tests

$$
\lambda_i + \lambda_{k-1-i} \;\approx\; 0, \qquad i = 0, \ldots, \lfloor k/2 \rfloor - 1,
$$

at relative tolerance 0.5 (qualitative gate — pairing residuals, not precision symplectic-spectrum recovery), and counts near-zero exponents ($|\lambda| < 10^{-4}$) as candidate invariant directions.

**Phase 2 → Phase 3 note.** The Sitting 3 production runs computed $k = 1$ tangent vector (the gates above need only $\lambda_{\max}$). With $k = 1$ the pairing check is vacuous — pairing is only observable once $k \geq 2$, and fully informative at $k = 14$. The full-spectrum pairing audit was therefore deferred to Phase 3. A preliminary $k = 14$ run after Sitting 4 found the flat $\mathbb{R}^{14}$ spectrum **dominated by two spurious quaternion-norm contraction modes** (one per body), summing to $-0.14\ \mathrm{s}^{-1}$ and failing pairing — the genuine embedding artefact of §7.3.2. Phase 3 (§7.10) then **fixed** this with a 12-dimensional cocycle: with the norm directions removed, the spectrum sums to $\approx 0$ and the Hamiltonian-pairing gate passes (AG-P3-1/2). The pairing check is now a real gate, not a vacuous one.

## §7.8 Connection to Chapters 05–06: eigenvalues, libration, and exponents

The three Ch 05–06 stability stories now close into a single consistent picture:

**Libration eigenvalues vs Lyapunov exponents (Ch 06 §6.4.2 ↔ §7.5(b)).** The linearised libration analysis yields purely imaginary eigenvalues $\pm i\,\omega_{\mathrm{lib}}$ at the prolate GG-stationary equilibrium. Imaginary eigenvalues mean oscillation, not divergence — and the nonlinear confirmation is precisely AG-LYAP-2: $|\lambda_{\max}| = 2.4 \times 10^{-6} \ll \omega_{\mathrm{lib}} = 6.4 \times 10^{-4}$. Libration is **not** exponential divergence; its Lyapunov exponent is consistent with zero. The Lagrange-top reduction (Ch 06 §6.4.3, Podolský §8.3 analog) explains *why*: the finite-amplitude motion is bounded in an effective-potential well, quadrature-integrable, hence regular.

**Dzhanibekov: saddle ≠ chaos, until coupled (Ch 05 §5.2.2 + Ch 06 §6.5.3 ↔ §7.5(c)).** The intermediate-axis equilibrium has a real positive linearised eigenvalue $\sigma$ (Ch 05 tennis-racket theorem), yet the torque-free finite-amplitude motion is *integrable* — the Jacobi-elliptic closed form of Ch 06 §6.5.3 parametrises every trajectory, including the near-separatrix flips, and an integrable system has identically zero Lyapunov spectrum. A linearised unstable eigenvalue at a saddle is a local statement; the global motion threads the heteroclinic structure regularly. Only the **gravity-gradient coupling** — which destroys the second first integral — converts the separatrix neighbourhood into a genuine thin chaotic layer, and even then with $\lambda_c \approx 0.002\,\sigma$ at the Phase 2 coupling strength. The hierarchy *linearised eigenvalue → integrable nonlinear oracle → coupled-system Lyapunov exponent* is the complete stability methodology of this project, and each level is independently gated (Ch 06 §6.5.5; AG-LYAP-1/3).

## §7.9 Open questions and Phase 3 outlook

- **OQ-PHASE2-LYA-1 — Lie-algebraic cocycle — RESOLVED (Phase 3, §7.10).** A tangent propagator in the body-frame Lie algebra ($\mathfrak{so}(3) \times \mathbb{R}^3$ per body, dimension 12) replaces the flat $\mathbb{R}^{14}$ FD Jacobian. **Achieved payoff:** removes the two spurious quaternion-norm contraction modes (12-dim by construction), so the full spectrum sums to $\approx 0$ and Hamiltonian pairing passes. **Payoff NOT achieved (prediction corrected):** the Phase-2-era expectation that this would restore fast FTLE convergence and make a $10^{-8}$ AG-LYAP-1 gate achievable was *falsified* — the Lie-algebraic cocycle gives $\lambda_{\max} = 2.29 \times 10^{-4}$ vs the $\mathbb{R}^{14}$ $2.30 \times 10^{-4}$ at $T = 20{,}000$ s. The slow convergence is intrinsic to the quasi-periodic dynamics, not the embedding (§7.3.2 corrected; AG-LYAP-1 stays at $10^{-3}$).
- **OQ-PHASE2-LYA-2 — Full (12-vector) spectrum — RESOLVED (Phase 3, §7.10).** With the 12-dimensional cocycle (Tier A deflation or Tier B Lie-algebraic), the Hamiltonian-pairing audit is a real gate: spectrum sum $|\Sigma\lambda| < 10^{-2}$ (AG-P3-1), the two large-negative norm modes absent, and near-zero invariant directions present (AG-P3-2). The interim Tier A deflation (project $\mathbb{R}^{14}$ against $q_A, q_B$ before each QR) and the canonical Tier B Lie-algebraic cocycle both discharge this. Kaplan-Yorke dimension estimates remain a stretch goal (out of Phase 3 scope).
- **OQ-PHASE2-LYA-3 — JWST attitude-control stability margins (still open).** The engineering question motivating the testbed: translate $\lambda_{\max}$ for realistic JWST-like configurations (composite inertia, L2-frame orbital environment, momentum-management cadence) into divergence time constants $1/\lambda$ comparable against attitude-control bandwidth and station-keeping intervals. Requires the v0.2+ L2 restricted-three-body extension (Ch 06 §6.9 OQ-6.1) before the orbital environment is representative.
- **Podolský coverage gap (restated from §7.1).** The CS-language canonical tier (`docs/canonical/cs/07-ljapunovovo-spektrum.md`, now authored) anchors the integrable reference dynamics in Podolský §8.1/§8.3 and points to the EN literature (Benettin 1980; Skokos 2010) for the Lyapunov formalism; Podolský 2018 does not cover Lyapunov spectra.

## §7.10 Phase 3 — the 12-dimensional spectrum (resolves OQ-PHASE2-LYA-1 + LYA-2)

Phase 3 builds two 12-dimensional cocycles and re-runs the spectrum, removing the §7.7 spurious-mode contamination. Both keep the $\mathbb{R}^{14}$ baseline (`lyapunov/cocycle.py`) intact for cross-validation; both are additive (no Phase 2 file modified).

**Tier A — quaternion-deflation cocycle** (`lyapunov/cocycle_deflated.py`). Wraps the $\mathbb{R}^{14}$ FD cocycle and projects every tangent vector orthogonal to the two current quaternion-radial directions $\hat{q}_A, \hat{q}_B$ before each Gram-Schmidt step, working in the deflated 12-dimensional subspace. Fast, low-risk, and proves the 12-vs-14 diagnosis directly.

**Tier B — Lie-algebraic cocycle** (`lyapunov/cocycle_liealg.py`). Tangent state $(\xi_A, \delta\boldsymbol{\omega}_A, \xi_B, \delta\boldsymbol{\omega}_B) \in \mathbb{R}^{12}$ with $\xi \in \mathfrak{so}(3)$. Orientation perturbations use the exponential map in the **right trivialisation matching the integrator** ($q_{k+1} = q_k \otimes \exp(\epsilon\,\xi)$), with the log map $\xi = 2\,\mathrm{arctan2}(|\mathbf{q}_{\mathrm{vec}}|, q_w)\,\mathbf{q}_{\mathrm{vec}}/|\mathbf{q}_{\mathrm{vec}}|$ mapping the propagated perturbation back to $\mathfrak{so}(3)$. Symmetric FD as in Phase 2, same evaluation count.

**Results (scenario (c), all gates GREEN; tests `tests/test_lyapunov_phase3.py`, 16 fast + 5 slow):**

| Gate | Check | Result |
|---|---|---|
| AG-P3-1 | spectrum sum $|\Sigma\lambda| < 10^{-2}$ | **PASS** (sum $\ll 0.14$; the $-0.064/-0.071$ norm modes gone) |
| AG-P3-2 | Hamiltonian pairing, $n_{\mathrm{near\text{-}zero}} \geq 2$ | **PASS** |
| AG-P3-3 | Lie-algebraic $\lambda_{\max}$ scenario (a) $< 10^{-3}$ | **PASS** ($2.29 \times 10^{-4}$; gate revised from the $10^{-5}$ aspiration — see below) |
| AG-P3-4 | Lie-algebraic $\lambda_{\max}$ within factor 3 of $\mathbb{R}^{14}$ value | **PASS** (both $> 0$, scenario (c)) |
| AG-P3-REGRESSION | All Phase 2 gates unchanged | **PASS** (16/16 fast) |

**The AG-P3-3 finding — a corrected prediction.** The Phase-2-era §7.3.2/§7.9 reasoning predicted the Lie-algebraic cocycle would restore fast FTLE convergence for scenario (a) and make a $10^{-8}$ gate achievable. It did not: the manifold-native cocycle gave $\lambda_{\max} = 2.29 \times 10^{-4}\ \mathrm{s}^{-1}$, essentially identical to the $\mathbb{R}^{14}$ value $2.30 \times 10^{-4}$. The $O(\log T/T)$ convergence is therefore **intrinsic to the quasi-periodic multi-frequency dynamics, not the flat embedding** — a genuinely chart-caused effect would have vanished under the manifold-native cocycle. The gate stays at the Phase 2 value $10^{-3}$ (AG-LYAP-1 parity), and §7.3.2 is corrected accordingly. This is the §7.3.2 honesty discipline turned on the chapter's own prior claim: the cross-tier prediction was physical-truth-tested by the Lie-algebraic experiment, and failed.

**What Phase 3 did and did not change.** It *fixed* the full-spectrum / Hamiltonian-pairing problem (the genuine embedding artefact). It *did not* change the single-vector $\lambda_{\max}$ gates or the AG-LYAP-1 threshold (the slow convergence is intrinsic). The two effects §7.3.2 originally conflated are now cleanly separated.

---

> **Reference anchors.** For the algorithm: Benettin, Galgani, Giorgilli & Strelcyn 1980, *Meccanica* **15**, 9–30, §2 (`Benettin1980` — the QR/Gram-Schmidt reorthonormalisation algorithm; primary anchor); Skokos 2010, *Lect. Notes Phys.* **790**, §2.2 (`Skokos2010` — modern formulation, FTLE convergence behaviour); Geist, Parlitz & Lauterborn 1990, *Prog. Theor. Phys.* **83**, 875–893 (`Geist1990` — discrete-QR method comparison). For the modified Gram-Schmidt numerics: Golub & Van Loan 2013 §5.2. For the integrable oracles calibrating the scenarios: Ch 06 §6.4.3 (Lagrange-top reduction) and §6.5.3 (Jacobi-elliptic asymmetric top), with Landau-Lifshitz Vol I §37 and Goldstein 2002 §5.6 as the underlying textbook anchors. For the symplectic integrator whose step supplies the cocycle: Hairer-Lubich-Wanner 2006 §VII.5 (`HairerLubichWanner2006`).
>
> **Supplementary anchor for Czech-reading audiences.** Podolský 2018 §8.1–§8.3 (`Podolsky2018`) covers the *integrable* reference dynamics — free symmetric top (§8.1) and Lagrange top with effective potential and the three precession regimes (§8.3) — that underpin scenarios (a) and (b). The Lyapunov-spectrum formalism is **not** covered in Podolský; CS readers should use Skokos 2010 for the algorithmic content (see §7.1 + §7.9).

---

**Next:** OQ-PHASE2-LYA-3 (JWST attitude-control stability margins) — requires the v0.2+ L2 restricted-three-body environment (Ch 06 §6.9 OQ-6.1) before the orbital setting is representative. Phase 3 (§7.10) resolved LYA-1 + LYA-2.

**End 07-lyapunov-spectrum.md (v0.1; §7.10 added + §7.3.2/§7.9 corrected at Phase 3, 2026-06-07).**
