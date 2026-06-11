# GLM memory-timescale pipeline

A hierarchy of increasingly mechanistic analyses of how recent **social motion**
shapes the moment-to-moment hazard of un-freezing. The logic is to *measure* the
memory timescale by model comparison (cheap, identifiable) and treat the temporal
kernel and the DDM as confirmation — not as the primary inference.

Every step runs on **one shared design** (`build_hazard_design.m`): the same
person-period intervals, events, baseline hazard, controls, and **bout-grouped CV
folds**, so all ΔCV-LL are paired and comparable. Each memory model is a single
z-scored predictor, so comparisons are multicollinearity-free and scale-free
(logistic CV-LL is invariant to linear rescaling of a predictor → only the feature
*shape* matters).

## Engine
- `build_hazard_design.m` — bouts + motion_cache → person-period design (records
  each interval's absolute frame so any causal filter can be sampled at it),
  baseline basis, controls, CV folds.
- `memory_feature.m` — one causal feature: `inst` | `box` (boxcar mean) |
  `exp` (EWMA = leaky-integrator state, leak τ) | `slope` (causal derivative).
- `fit_cv.m` — logistic hazard fit → AIC + grouped CV-LL (+per-fold) + in-sample LL.
- `profile_timescale.m` — CV-LL vs timescale grid for `box`/`exp`.
- helpers: `geto.m`, `raised_cosine_basis.m`, `grouped_cv.m`.

## Steps
1–5 are in **`memory_timescale_analysis.m`** (run this first):
1. **motion matters** — baseline vs baseline+mean (ΔCV-LL, ΔAIC).
2. **temporal structure** — mean vs mean+slope (derivative).
3. **boxcar timescale** — CV-LL(T), the averaging horizon.
4. **exponential timescale** — CV-LL(τ); best τ = the **leak time-constant**.
5. **memory architecture** — instantaneous vs best boxcar vs best exp.

Queued (separate scripts, reuse the engine + same folds):
6. **kernel as confirmation** — smoothed `run_hazard_kernel`, overlay best-exp / best-box implied kernels.
7. **bootstrap identifiability** — cluster bootstrap over bouts → distribution of best-T / best-τ / win margins.
8. **shuffle nulls** — circular-shift sm (keeps autocorrelation, breaks alignment): does the timescale survive?
9. **mechanistic DDM** — generative leaky-accumulator fit to durations + sm (done last).

## Notes
- `sm` is globally z-scored; missing → baseline (0); filters run on the *full*
  fly trace so memory has equilibrated by freeze onset and integrates pre-freeze history.
- Right-censored bouts (collision) contribute at-risk intervals but no event.
- Results/figures save to `paths.dataset` as `glm_memory_timescale_<stamp>.{mat,png}`.
