# recovery — accumulator-driven datasets → hazard-kernel recovery

Parameter-recovery validation for the social-motion hazard-kernel analysis.

We generate **synthetic freeze datasets** in which the freeze duration is the
first-passage time of a **leaky accumulator** integrating the *real* social-motion
each fly experienced, then run the **same** hazard-kernel analysis used on real
data (`run_hazard_kernel`, the function-form core of
`../hazard_kernel_sm.m`) and check whether it recovers the leak
time-constant we put in.

Because the accumulator is driven by the real `sm` trace — the same series the
analysis later reads back — the regressor autocorrelation matches the real
analysis, so this is an honest test of the method rather than a toy.

## Files

- **`generate_accumulator_dataset.m`** — given the real bout template + per-bout
  `sm` drive + generative params `P`, simulates `P.n_aug` accumulator
  trajectories per real bout (augmentation) and returns a synthetic bout table
  (`bl`) ready for `run_hazard_kernel`, plus duration/censoring diagnostics.
  Uses `sim_leaky_accumulator` (`../../../sims/simulators`).
- **`kernel_recovery_sweep.m`** — driver. Loads real bouts + `motion_cache`,
  builds the `sm` drive, and sweeps one generative parameter (`sweep_var`):
  - `'beta'` — social GAIN, leak fixed at `tau_fixed` (specificity / amplitude
    test). Plots recovered kernel amplitude vs `beta`; should grow with `beta`
    and collapse to ≈0 at `beta = 0` (sm no longer drives the accumulator).
  - `'tau'` — leak TIME-CONSTANT, gain fixed (timescale recovery). Plots
    recovered `tau` vs generative `tau` against the identity line.
  (Optionally calibrates `theta` to the real median.) Saves everything
  (synthetic datasets + recovered kernels) to `results/`.
- **`kernel_temporal_structure.m`** — a *different* question: does the temporal
  structure of `sm` matter? Holds the generative params fixed and varies only
  how the drive is **constructed** — full time-resolved `sm(t)` vs the bout-mean
  scalar (`mu0 + β·mean(sm_bout)`, constant drift), with optional `white`/`ar1`
  drives. Headline readout is the **acausal/causal kernel-weight ratio**: ≈0 ⇒
  causal/time-resolved, ≈1 ⇒ flat & symmetric ⇒ only the level matters (timing
  irrelevant). The `mean` condition still lets the analysis read the real
  time-resolved `sm`, so it shows what the kernel does when no timing was used.

## Run

```matlab
kernel_recovery_sweep        % edit the params block at the top first
```

## The leak ↔ kernel mapping

| generative leak | `P.lambda` | expected recovered kernel |
|---|---|---|
| `tau = Inf` | `0` | flat (perfect accumulator, infinite memory) |
| `tau` finite | `1/tau` | exponential, timescale ≈ `tau` (leaky integrator) |
| `tau → dt` | `1/dt` | delta at `τ=0` (extrema / instantaneous) |

## Knobs (`P` and `opts` in the driver)

- `sweep_var` (`'beta'` or `'tau'`) + `sweep_list` choose what is swept;
  `tau_fixed` is the leak used during a `beta` sweep, `P.beta` the gain used
  during a `tau` sweep.
- `P.mu0` baseline drive (sets the duration timescale), `P.sigma` diffusion SD,
  `P.n_aug` trajectories per real bout, `P.theta` bound.
- `calibrate` matches `theta` to the real median per swept value (keeps
  durations realistic).
- `opts` mirrors `hazard_kernel_sm.m` (kernel window, basis, penalty). Set
  `opts.n_shuffle > 0` and `opts.do_stage2 = true` to also run the temporal-order
  shuffle control and the integration-vs-extrema model comparison (slower).

## Reading the result

- **Fig 2** is the headline, and depends on `sweep_var`:
  - `beta` sweep — recovered kernel **amplitude** (peak and Σ of the causal
    kernel) vs generative `beta`: should rise with `beta` and sit at ≈0 for
    `beta = 0` (the null / specificity check).
  - `tau` sweep — recovered `tau` vs generative `tau` against the identity
    line; `τ_com` (model-free center-of-mass) and `τ_exp` (exponential fit)
    should both track.
- **Fig 1** overlays the recovered kernels — they should broaden as `tau` grows
  (sharp/peaked for short `tau`, flat for the accumulator).
- **Fig 3** checks the synthetic durations resemble the real ones (KS in console).
- Console also prints censoring %, the trunc-filter loss, and `R²` of the
  exponential leak fit per `tau`.

## Caveats

- Augmented replicates of one real bout share an `sm` history; they are grouped
  under a shared `bout_id` so grouped CV does not leak them across folds.
- `theta` calibration matches the *median* only — the full shape match is the KS
  diagnostic, not a fit target.
- The recovery readouts (`τ_com`, `τ_exp`) come from the kernel **shape** and do
  not depend on CV; `opts.n_shuffle`/`do_stage2` only affect the optional
  temporal-order and model-comparison outputs.
