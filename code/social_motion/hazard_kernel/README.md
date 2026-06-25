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
