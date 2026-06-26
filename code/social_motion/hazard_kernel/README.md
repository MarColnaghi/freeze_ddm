# hazard_kernel — social motion → freeze-termination hazard

All of this folder asks one question: **does the social motion a focal fly
experiences drive *when* it un-freezes, and over what *memory timescale*?**
Three methods attack it on real data, plus a synthetic parameter-recovery layer
that validates the kernel method.

## Shared core
- **`run_hazard_kernel.m`** — discrete-time hazard kernel of freeze termination
  vs social motion. The analysis core (extracted from `hazard_kernel_sm.m`) so the
  real data and the synthetic ground truth run through *identical* methodology.
  Used directly by everything in `recovery/`.

## Methods (on real data)
- **`hazard_kernel_sm.m`** — the time-resolved kernel driver:
  `logit P(unfreeze│frozen) = baseline(t) + Σ_τ β(τ)·sm(t−τ) + controls`.
  The kernel *shape* distinguishes integration (broad/sustained) from
  instantaneous (sharp near τ≈0); a small acausal segment is a reverse-causation
  control. *(Note: currently re-implements `run_hazard_kernel`'s logic inline
  rather than calling it — keep them in sync.)*
- **`hazard_model.m`** — a leaner discrete-time hazard (0.1 s bins, 30 lags,
  leaky integrator). Same idea, fewer moving parts.
- **`glm_memory/`** — a model-comparison pipeline that *measures* the memory
  timescale cheaply/identifiably (boxcar vs EWMA features, grouped CV-LL), with
  the kernel and DDM as confirmation. See `glm_memory/README.md`.

## Validation (on synthetic data)
- **`recovery/`** — generates synthetic freezes whose durations are leaky-
  accumulator first-passage times driven by the *real* `sm`, then runs
  `run_hazard_kernel` on them and checks it recovers the leak time-constant /
  gain we put in. See `recovery/README.md`.

## Layout
```
hazard_kernel/
├── run_hazard_kernel.m        shared estimator (called by recovery/)
├── hazard_kernel_sm.m         real-data kernel driver
├── hazard_model.m             leaner hazard variant
├── glm_memory/                model-comparison memory-timescale pipeline
└── recovery/                  synthetic parameter-recovery (+ results/)
```

## Relation to reverse-correlation kernels (Okazawa et al. 2018)
Okazawa, Sha, Purcell & Kiani 2018 (*Nat. Commun.* 9:3479) also speak of a
"kernel," but it is a **different object** from our `β(τ)`:

- **Their psychophysical kernel** is a *reverse correlation*:
  `K(t) = E[s(t)|choice1] − E[s(t)|choice2]` — average the stimulus *conditioned
  on the outcome*, subtract the two choices. It targets the linear **sensory
  weight** `w(t)` (`K ∝ w` under SDT/simple DDM), and their whole point is that
  it is *distorted* by the decision process (bound crossing, non-decision time).
  One binary choice per trial; x-axis = position within trial / time-to-response.
- **Our hazard kernel** `β(τ)` is a *forward* penalized-logistic GLM on
  person-period survival data: `logit P(unfreeze|frozen) = baseline(t) +
  Σ_τ β(τ)·sm(t−τ) + controls`. It targets the **integration / memory timescale
  itself** (broad ⇒ integration, sharp near τ≈0 ⇒ instantaneous). Many at-risk
  steps per bout, one absorbing event; x-axis = lag τ relative to "now".

So the two are duals: `E[stim|outcome]` (reverse) vs. `P(event|stim history)`
(forward). What Okazawa treat as a *nuisance* (decision dynamics) is exactly our
*target*; their DDM-simulation validation is mirrored by `recovery/`.

**Bridge — STA↔GLM duality.** For white/Gaussian input an event-triggered
average and a forward-GLM filter coincide up to scale. We already compute *both*:
`eta_exc` (`run_hazard_kernel.m`) is an Okazawa-style reverse-correlation kernel
on our data, and `rho_eta` correlates it against the fitted `β(τ)`. They diverge
when the input is autocorrelated — which `sm` is — and only the forward GLM
partials that out. This maps Okazawa's "kernel ≠ sensory weight" onto our
"event-triggered sm ≠ causal hazard kernel"; our acausal τ<0 segment, time-shuffle
null, and ACF diagnostic (`rho_acf`) are the freeze-DDM analog of their
non-decision-time correction.

**Practical note.** Their non-decision-time result predicts our *causal* kernel
may peak at a short *positive* lag (sensory/motor latency), not exactly at τ=0 —
a dip at τ≈0 is expected, not evidence against integration.
