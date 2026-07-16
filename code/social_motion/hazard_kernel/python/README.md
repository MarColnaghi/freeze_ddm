# Social-motion hazard kernel — Python port

A compact, self-contained NumPy/Matplotlib reimplementation of the MATLAB
hazard-kernel pipeline (`../visual_aid/kernel_fit_single.m`,
`../visual_aid/kernel_recovery_ddm.m`).

The model is a discrete-time survival / person-period logistic GLM for the
un-freeze hazard, driven by a temporal kernel over recent social motion:

```
logit P(break at interval i | still frozen) = b0(time-in-freeze) + sum_tau beta(tau) * sm(t - tau)
```

`beta(tau)` is expanded in a small raised-cosine basis (`beta = B @ b`), fit by
penalised IRLS (Newton) with an optional 2nd-difference smoothness penalty
`lambda * D'D`, and `lambda` is chosen by bout-grouped cross-validation.

## Files

- `hazard_kernel.py` — the library (numpy only):
  - `raised_cosine_basis`, `second_difference` — basis `B` and penalty operator `D`.
  - `penalized_logit(X, y, P)` — IRLS logistic with penalty `P`, returns `(b, V)`.
  - `grouped_cv_loglik`, `select_lambda` — bout-grouped CV to pick `lambda`.
  - `simulate_accumulator` — leaky-accumulator (DDM) first-passage simulator
    (`lam = 0` → perfect DDM / flat kernel; `lam = 1/tau` → exponential kernel).
  - `build_design` — stack every at-risk interval into the long design.
  - `fit_leak_tau` — coarse exponential fit to read off the recovered timescale.
- `recover_demo.py` — end-to-end **ground-truth recovery**: simulate un-freeze
  times from a DDM with a KNOWN kernel, recover `beta(tau)`, and overlay the two.
  Saves `kernel_recovery_ddm_py.png` (CV curve + recovered-vs-truth kernel).
- `recover_sklearn.py` — same recovery, but fit with scikit-learn's **built-in
  regularizers** (`LogisticRegression`, `l1_ratio` 0 = ridge / 1 = lasso), with
  strength `C` chosen by bout-grouped CV (`GroupKFold` + `GridSearchCV`,
  `scoring='neg_log_loss'`, `StandardScaler` in a pipeline). Overlays sklearn
  ridge (L2), sklearn lasso (L1), the custom curvature P-spline, and ground
  truth. Saves `kernel_recovery_sklearn_py.png`.

## Run

```bash
python recover_demo.py
```

Knobs at the top of `recover_demo.py`: `TAU_GEN` (`np.inf` = flat/perfect DDM),
`BETA`/`SIGMA`/`THETA` (generative drive/noise/bound), `N_BOUTS`/`N_AUG`
(number of events), `NB_BASIS`, `CV_FOLDS`.

## What to check (the point of the recovery)

- **Shape:** the recovered kernel tracks the exponential ground truth on the
  causal side.
- **Timescale:** printed `tau_hat` is in the right ballpark of `tau_gen`
  (biased low — the bound + diffusion noise make the tau→kernel mapping only
  approximate; this is expected).
- **Acausal:** the future side stays near 0 (a small blur just past `tau=0` is
  the smooth basis being unable to represent the sharp causal cutoff).
- **CV curve:** nearly flat across `lambda` then drops — smoothness barely
  changes prediction, so `lambda` is only weakly identified (same finding as on
  the real data).

## Environment note

The repo's base Anaconda has a broken NumPy (2.x) / SciPy·Matplotlib (built for
1.x) ABI mismatch, so `import matplotlib` fails there. Use an isolated env:

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python recover_demo.py
```

`recover_demo.py` needs only `numpy` + `matplotlib` (the penalised IRLS and CV
are implemented directly, mirroring the MATLAB). `recover_sklearn.py` also needs
`scikit-learn`.

## Built-in regularizers vs the curvature penalty

`recover_sklearn.py` uses scikit-learn's built-in penalties. Worth knowing what
they can and cannot express:

- **Ridge (L2, `l1_ratio=0`)** penalises `sum(b^2)` — shrinks bump coefficients
  toward 0 (magnitude).
- **Lasso (L1, `l1_ratio=1`)** penalises `sum|b|` — drives some bumps to exactly
  0 (sparse).
- **P-spline (custom IRLS)** penalises `||D b||^2` — the *curvature* of the
  coefficient sequence, i.e. smoothness, without shrinking the overall level.

sklearn's built-ins act on coefficient **magnitude**, applied **uniformly** to
every coefficient — they cannot express the curvature penalty, nor exempt the
baseline block. That is exactly why the production pipeline uses a hand-built
IRLS with a block-specific smoothness penalty. In the strong-signal recovery
demo all three land on essentially the same kernel (the 12-bump basis already
supplies most of the smoothing); the regularizer choice matters more when events
are scarce and the fit is noisier.
