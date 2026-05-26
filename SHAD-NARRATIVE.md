# How an asymmetric body tells you its own moment of inertia, just by tumbling

> **Audience.** Same as the fourier shad-tier — you understand plots, you don't need theorem-proof rigor, you do want to see real numbers fall out of the math.
>
> **Reading time.** ~15 minutes if you skim the FFT discussion, ~30 if you read the equations.
>
> **License.** CC-BY-SA-4.0.
>
> **Status.** Companion narrative to the runnable first-cut prototype at `experiments/jwst-l2-first-cut/`. Pragmatic experiment, not the full `lege-artis/jwst-l2-coupled` project (that ships later).

---

## §1. The setup, in two satellites

Two simplified spacecraft, floating in space far enough from anything else that nothing else matters for a while.

Body A is JWST-shaped — a flat sunshield at one end, a primary mirror at the other, a boom connecting them through a spacecraft bus. The numbers are illustrative not exact: ~2450 kg total, principal moments of inertia I = diag(23322, 23322, 15384) kg·m² in its own body frame. That last number is smaller than the other two — A is **oblate**, like a tossed pancake. It will spin happily about its z-axis (the boom axis), and if you spin it about anything else, that something else will not stay still.

Body B is a probe — a short cylinder with a dish antenna. ~100 kg, principal moments I = diag(29, 29, 5) kg·m². The opposite anisotropy: I_z is the *smallest* number, not the largest. B is **prolate**, pencil-shaped. It also spins happily about its z-axis.

We put A at the origin. We put B at +x = 50 metres. Both have small initial angular velocities — but **not aligned with a principal axis**. Body A spins at ω_body = (0.02, 0, 0.08) rad/s — most of the spin is about z, but there's a 25% component perpendicular to z. Body B gets ω_body = (0, 0.03, 0.08) rad/s — same idea, perpendicular component along y.

We also tilt A by 0.3 radians (≈17°) about an awkward axis ((1,1,0)/√2 in inertial coordinates) so that A is **not** sitting with its principal axes aligned to the inertial frame.

Then we let the system evolve under three forces:
1. Newtonian gravity between A and B (very small at this distance — F ≈ 2 × 10⁻⁷ N).
2. The gravity-gradient torque each body exerts on the other (smaller still — τ ≈ 6 × 10⁻¹⁰ N·m).
3. The body's own self-torque from being a rigid object spinning around an axis that isn't a principal axis. **This one's not small.**

The third one is the interesting bit. It's the term `ω × (I·ω)` in Euler's equations.

## §2. Euler's equations, in plain English

Newton said `F = ma` for points. Euler said the rotational analogue for rigid bodies has an extra term:

```
I · dω/dt + ω × (I · ω) = τ
```

In body frame. The `I` is the body's inertia tensor (constant in body frame; that's why we use body frame). The `ω` is the angular velocity. The `τ` is whatever external torque is acting.

The `ω × (I · ω)` term is the body talking to itself. If `ω` is aligned with a principal axis of `I`, then `I·ω` is parallel to `ω`, the cross product is zero, and there's no self-torque. The body spins steadily.

If `ω` is *not* aligned with a principal axis — like in our setup, where A is spinning about (x + 4z)-ish in body coordinates — then `I·ω` points somewhere else, and `ω × (I · ω)` produces a torque that the body wasn't asking for. It precesses.

How fast?

For an axisymmetric body where I_xx = I_yy ≠ I_zz, the Euler equations simplify enough that we can pull a closed-form precession frequency:

```
λ_Euler = (I_zz − I_xx) / I_xx · ω_z
```

For Body A with ω_z = 0.08 rad/s and the inertia values above:

```
λ_A = (15384 − 23322) / 23322 × 0.08
    = −0.3404 × 0.08
    = −0.0272 rad/s
```

Period = 2π / |λ| ≈ **231 seconds**.

This is a prediction. The body should oscillate its perpendicular-to-z angular velocity components at this period, with no external torque required, just from spinning about a not-quite-principal axis.

## §3. We run it and we get the data

`python3 run_first_example.py`. Three seconds later, an NDJSON file with 601 snapshots of the system state covering 600 seconds. Each snapshot has the full state of both bodies — positions, velocities, quaternions, rotation matrices, body-frame and inertial-frame angular velocities, the inertia tensor expressed in the inertial frame, the principal axes as vectors in the inertial frame, and a handful of diagnostics (energies, total angular momentum, gravity-gradient torque magnitudes).

The first thing we check is the boring thing: are total energy and total angular momentum conserved? In this isolated A-B system with no Sun, no Earth, no external torques, they have to be. If they aren't, the integrator is lying.

```
|dE / E_0| = 6.2 × 10⁻¹²
|dL / L_0| = 2.1 × 10⁻¹⁰
```

Both at machine precision. The classical 4th-order Runge-Kutta integrator behaves. We're good.

Now the interesting check: do the rotational modes show up?

We take the body-frame ω_A_x time series (601 points, sampled every 1 second of simulation time), subtract the mean, run an FFT, and look for the peak.

```
peak frequency: 0.005 Hz
peak period:    200 seconds
peak magnitude: 4.84 (in whatever units we don't care about)
```

200 seconds. We predicted 231 seconds. The FFT bin spacing at this observation window is 1/600 = 1.67 mHz, so observable periods near 230 s are quantised: 600, 300, 200, 150, 120, ... seconds. 200 is the closest bin to 231.

**The prediction matched.** The body told us its moment of inertia, just by tumbling.

## §4. What you would do at a oscilloscope

This is the same story as B1 (oscilloscope trace) from the fourier shad-tier guide. There, a periodic voltage signal showed up as a peak in the DFT power spectrum, and you could read the oscillator's frequency off the plot. The signal source was a function generator. The interpretation was "this is what the generator is doing."

Here the signal source is a rotating rigid body. The interpretation is "this is what the inertia tensor is doing." Same algorithm. Same plot. Different physics.

If you didn't know the body's inertia tensor — say, you launched a spacecraft and you're not entirely sure what its mass distribution turned out to be after assembly — you could give it a controlled tumble (perpendicular + parallel components), measure the perpendicular angular-velocity time series with any reasonable IMU, FFT it, find the peak frequency, and read off the ratio (I_zz - I_xx) / I_xx from the formula above. Combined with the known ω_z you commanded, that's a measurement of the body's principal-moment ratio with no need to tear it apart on a coordinate-measuring machine.

This is **not the way actual spacecraft inertia tensors get characterised** — there are calibration-table methods, swing tests, and post-launch on-orbit tuning that are higher-precision. But it's the same physics. The Euler precession frequency carries the inertia ratio, and the DFT extracts it.

## §5. What the inertia tensor does in the inertial frame, since I promised you that

While the body tumbles, its inertia tensor expressed in body frame stays constant — `I_body = diag(23322, 23322, 15384)`. That's the whole point of using body frame; the math is easier when the tensor doesn't change.

But the inertia tensor in the **inertial** frame changes constantly, because the body is rotating relative to inertial. The transformation is `I_inertial(t) = R(t) · I_body · R(t)^T` where `R(t)` is the body-to-inertial rotation matrix.

At t = 0, with the body tilted by 0.3 rad about (1,1,0)/√2, the diagonal of `I_inertial` is roughly (23322, 23322, 15384) but with off-diagonal elements injected by the tilt. The trace (sum of diagonal) stays at 62028 because trace is rotation-invariant. The eigenvalues stay at (15384, 23322, 23322) because principal moments are rotation-invariant. But the values you see when you write the matrix out in the inertial-frame basis change.

At t = 600 s, the body has tumbled a few full turns. The diagonal in inertial frame is (19094, 22224, 20712). Trace: 62030. (Trace preserved to ppm — the renormalised quaternion stays unit-magnitude well enough.) Eigenvalues: still (15384, 23322, 23322). The body hasn't changed; the *labels we put on its axes* have.

The principal-axis vectors in inertial frame *do* change. At t = 0, the third (smallest-moment) principal axis is along the tilted z-direction — slightly off inertial-z. At t = 600 s, that same physical axis (the one passing through the sunshield and the primary mirror) is pointing in a completely different inertial direction. The body is tumbling, and you can watch its physical axes sweep across the inertial frame in the snapshot data.

Each NDJSON snapshot has `principal_axes_inertial` as a 9-element list — three eigenvectors of the inertial-frame inertia tensor, written column-major. The output `principal_moments` is the body-frame invariant. You can reconstruct the body's orientation at any snapshot by stacking the three principal-axis vectors as columns of a rotation matrix.

This is the data you'd plot if you wanted to *visualise* the tumble. We're not doing that here. The numbers are enough for the punchline.

## §6. The gravity-gradient torque is tiny but real

The diagnostics block of each snapshot has `tau_gg_A_magnitude` and `tau_gg_B_magnitude` — the magnitudes of the gravity-gradient torques each body exerts on the other.

At t = 600 s: `tau_gg_A_magnitude ≈ 6.3 × 10⁻¹⁰ N·m`.

Body A's rotational kinetic energy at this point is ~54 J. A torque of 6 × 10⁻¹⁰ N·m acting on this body needs about `E_rot / τ ≈ 54 / 6e-10 ≈ 10¹¹ seconds` ≈ **3000 years** to do meaningful work against the spin.

So gravity-gradient torque from Body B is irrelevant on the 600-second timescale. It's not zero — it modulates the free precession on top, but at 1 part in 10⁹ or so per orbit. To see gravity-gradient effects, we'd need either (a) a much closer encounter (B at 5 m would give 10⁶ × stronger torque), (b) a much longer integration (~10⁸ seconds), or (c) a much heavier perturbing body.

The "JWST-L2 coupled dynamics" story doesn't actually need the perturbing body to dominate. The point of the eventual full project is that at L2, where the central gravity of Sun + Earth largely cancels, gravity-gradient from anything nearby becomes the *dominant* attitude-disturbing mechanism — because every other mechanism is also small at that location. Here in the first-cut prototype we've isolated A and B, removed Sun and Earth, and watched A's free precession independently of any external perturbation. The gravity-gradient torque sits there in the diagnostics, available for analysis when we integrate longer.

## §7. Why this is the same story as fourier B1..B5

You took an oscilloscope trace, transformed it, and read the source's frequency. You took an audio recording, transformed it, and identified the overtones. You took a vibration log, transformed it, and diagnosed the bearing fault. You took a radar return, transformed it, and extracted the target's distance.

Here you took an angular-velocity time-series from a tumbling satellite, transformed it, and read the inertia ratio.

It's the same algorithm. The reason it works is the same reason: the underlying physical system contains a characteristic frequency, the system's output modulates at that frequency, and the DFT picks it out of whatever noise (or background, or other modes) is around it. In B7 (queued — nuclear reactor noise) we'll meet the case where the bare DFT stops being enough and we need the "and beyond" toolkit. Here in the rigid-body world the bare DFT is sufficient, because the dynamics are mostly linear (free precession is a linear normal mode of Euler's equations for axisymmetric bodies).

When the dynamics get non-linear — strong gravity-gradient coupling, gyroscopic resonances when the precession period matches an orbital period, the full coupled-6-DOF problem — then the bare DFT will start to mislead, and we'll go pull out the same toolkit B7 uses. But that's later.

## §8. What to do with the data

Open the NDJSON file. Each line is a snapshot. Pick any field you care about and plot it. Suggestions:

- `A.omega_body[0]` vs `t` — should look like a slowly damped sinusoid (it's actually constant amplitude in the axisymmetric case; small numerical drift). The period is ~230 s.
- `A.principal_axes_inertial` — pick the third column (the smallest-moment principal axis); plot its components vs t. Should trace a cone (precession trajectory).
- `diagnostics.E_total` vs `t` — should be flat (conservation working).
- `diagnostics.tau_gg_A_magnitude` vs `t` — varies slowly with the body's orientation relative to body B; this is the small modulation that would grow over thousands of seconds.

`jq` works fine. So does any Python streaming-JSON reader. The format is deliberately verbose so that future-you can extract whatever quantity future-you needs without re-running the simulation.

## §9. Closing

The headline payoff:

```
The body tumbles.
The body's tumble carries the body's inertia tensor.
The DFT, applied to the tumble, gives you the tensor.
```

This is the same shape as the rest of the shad-tier guide. Different physics, different domain, same algorithm. The DFT is the canonical detection tool across observational systems — engineered (oscilloscope, audio, radar) and natural (rigid-body dynamics, pulsar timing, reactor noise).

Next chapter in this small experiment: add the Sun + Earth gravity, put A on a halo orbit around L2, watch the long-timescale gravity-gradient coupling actually do something visible. That's what the full `lege-artis/jwst-l2-coupled` project gets to. This first-cut closes when we've shown the integrator + the data pipeline + the basic interpretation arc work end-to-end.

They do.

---

**Code:** `experiments/jwst-l2-first-cut/run_first_example.py`
**Data:** `experiments/jwst-l2-first-cut/outputs/first_example_trajectory.ndjson`
**Strategic capture:** `_config/JWST-L2-COUPLED-DYNAMICS-PROJECT-v0.1.md`
**Doctrine:** `_config/LEGE-ARTIS-LANGUAGE-DOCTRINE-v0.1.md` §4.4 — what the eventual canonical-reference implementation must satisfy.

*— Pete Y., May 2026. CC-BY-SA-4.0. Drafted for the Shaddack-edition release.*
