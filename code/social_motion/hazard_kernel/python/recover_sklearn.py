"""
recover_sklearn.py — same hazard-kernel recovery, but with scikit-learn's
BUILT-IN regularizers, contrasted with the custom curvature (P-spline) penalty.

We reuse the simulator + design builder from hazard_kernel.py, then fit beta(tau)
four ways and overlay them on the DDM ground truth:

  * sklearn L2  (ridge)  — penalises sum(b^2): shrinks bump coefficients toward 0.
  * sklearn L1  (lasso)  — penalises sum|b|:  drives some bumps to exactly 0 (sparse).
  * P-spline  (custom)   — penalises ||D b||^2: curvature -> smooth kernel (our IRLS).
  * ground truth exp(-tau/tau_gen).

Regularisation strength is chosen by BOUT-GROUPED cross-validation
(GroupKFold + GridSearchCV, scoring = held-out log-loss), matching the MATLAB.

Key teaching point: sklearn's built-in penalties act on coefficient MAGNITUDE
(ridge/lasso), applied uniformly to every coefficient. They cannot express the
CURVATURE penalty (||D b||^2) or exempt the baseline block — which is exactly why
the real pipeline uses a hand-built IRLS with a block-specific smoothness penalty.

Run:  python recover_sklearn.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import hazard_kernel as hk


# ── config (mirror recover_demo.py) ──────────────────────────────────────────
FPS = 60
KERNEL_PAST_S, KERNEL_FUTURE_S = 4.0, 1.0
DT_FRAMES, ENTRY_FR = 6, 30
NB_BASIS, BUMP_WIDTH, NB_BASELINE = 12, 1.25, 8
CV_FOLDS, GRID_ANCHOR, MASK_PREONSET = 5, "offset", True

TAU_GEN = 1.0
MU0, BETA, SIGMA, THETA = 5.0, 0.0, 1, 7.0
DT, MAXT_FR = 1.0 / FPS, int(round(10.5 * FPS))
N_BOUTS, N_AUG, SEED = 2000, 5, 7


def make_social_motion(n_bouts, rng):
    out = []
    for _ in range(n_bouts):
        length = int(rng.integers(180, MAXT_FR))
        z = rng.standard_normal(length)
        a, x = 0.94, np.zeros(length)
        for k in range(1, length):
            x[k] = a * x[k - 1] + np.sqrt(1 - a * a) * z[k]
        out.append((np.abs(x) + 0.3 * rng.random()).astype(float))
    return out


def grouped_gridsearch(X, y, groups, l1_ratio):
    """Pick C by bout-grouped CV (held-out log-loss), StandardScaler + saga.

    l1_ratio = 0 -> pure L2 (ridge), l1_ratio = 1 -> pure L1 (lasso) — the current
    sklearn API (the old penalty='l2'/'l1' is deprecated). Coefficients are mapped
    back out of the StandardScaler so beta = B @ coef is in the original units.
    Returns (coef_original, C).
    """
    Cs = np.logspace(-4, 4, 15)
    gkf = GroupKFold(n_splits=CV_FOLDS)
    pipe = make_pipeline(
        StandardScaler(),
        LogisticRegression(l1_ratio=l1_ratio, solver="saga", max_iter=5000, tol=1e-3))
    gs = GridSearchCV(pipe, {"logisticregression__C": Cs},
                      cv=gkf, scoring="neg_log_loss", n_jobs=-1)
    gs.fit(X, y, groups=groups)
    best = gs.best_estimator_
    lr = best.named_steps["logisticregression"]
    sc = best.named_steps["standardscaler"]
    coef = lr.coef_.ravel() / sc.scale_          # undo standardisation -> original space
    return coef, gs.best_params_["logisticregression__C"]


def main():
    rng = np.random.default_rng(SEED)

    # 1. templates + augmentation --------------------------------------------
    templates = make_social_motion(N_BOUTS, rng)
    sm_list, bout_of_row = [], []
    for _ in range(N_AUG):
        for b in range(N_BOUTS):
            sm_list.append(templates[b]); bout_of_row.append(b)

    # 2. simulate DDM freezes -------------------------------------------------
    lam = 0.0 if np.isinf(TAU_GEN) else 1.0 / TAU_GEN
    drives = [MU0 + BETA * (sm / 10.0) for sm in sm_list]
    durations, censored = hk.simulate_accumulator(drives, THETA, SIGMA, lam, DT, seed=SEED + 100)

    # 3. lag grid + design ----------------------------------------------------
    back, fwd = int(round(KERNEL_PAST_S * FPS)), int(round(KERNEL_FUTURE_S * FPS))
    lag_fr = np.arange(-fwd, back + 1)
    lag_s = lag_fr / FPS
    caus = lag_fr >= 0
    allsm = np.concatenate(templates)
    sm_mu, sm_sd = allsm.mean(), allsm.std()

    S, y, tinf, gid_local = hk.build_design(
        durations, censored, sm_list, lag_fr, FPS, DT_FRAMES, ENTRY_FR,
        GRID_ANCHOR, MASK_PREONSET, sm_mu, sm_sd)
    gid = np.array([bout_of_row[i] for i in gid_local])
    print(f"Design: {S.shape[0]} intervals, {int(y.sum())} events (hazard {y.mean():.3f})")

    # 4. basis + baseline block ----------------------------------------------
    B = hk.raised_cosine_basis(lag_fr, NB_BASIS, BUMP_WIDTH)
    K = B.shape[1]
    Xk = S @ B
    ut = (np.argsort(np.argsort(tinf)) + 1) / tinf.size
    Bb = hk.raised_cosine_basis(ut * (back + fwd) - fwd, NB_BASELINE, 1.25)[:, 1:]
    Xsk = np.column_stack([Xk, Bb])           # sklearn adds its own (unpenalised) intercept
    kcols = np.arange(K)

    # 5a. sklearn built-in L2 (ridge) and L1 (lasso) --------------------------
    coef_l2, C_l2 = grouped_gridsearch(Xsk, y, gid, l1_ratio=0.0)   # ridge
    coef_l1, C_l1 = grouped_gridsearch(Xsk, y, gid, l1_ratio=1.0)   # lasso
    kernel_l2 = B @ coef_l2[kcols]
    kernel_l1 = B @ coef_l1[kcols]
    nnz = int(np.sum(np.abs(coef_l1[kcols]) > 1e-8))
    print(f"sklearn ridge (l2): C={C_l2:.3g}")
    print(f"sklearn lasso (l1): C={C_l1:.3g}, {nnz}/{K} kernel coefs nonzero")

    # 5b. custom P-spline (curvature) via our IRLS ----------------------------
    Xp = np.column_stack([np.ones(S.shape[0]), Xk, Bb])
    ker = np.arange(1, 1 + K)
    D = hk.second_difference(K)
    Pk = np.zeros((Xp.shape[1], Xp.shape[1])); Pk[np.ix_(ker, ker)] = D.T @ D
    lam_grid = np.logspace(-1, 5, 10)
    lam_p, _ = hk.select_lambda(Xp, y, Pk, gid, lam_grid, CV_FOLDS, seed=SEED, rule="argmax")
    bP, _ = hk.penalized_logit(Xp, y, lam_p * Pk)
    kernel_p = B @ bP[ker]
    print(f"P-spline (curvature): lambda={lam_p:.3g}")

    # 6. ground truth ---------------------------------------------------------
    gt = np.zeros_like(kernel_p)
    gt[caus] = 1.0 if np.isinf(TAU_GEN) else np.exp(-lag_s[caus] / TAU_GEN)

    def scale_to(k):  # least-squares scale ground truth onto a recovered kernel (causal)
        a = (gt[caus] @ k[caus]) / max(gt[caus] @ gt[caus], 1e-12)
        return a * gt

    for name, k in [("P-spline", kernel_p), ("ridge", kernel_l2), ("lasso", kernel_l1)]:
        r = np.corrcoef(k[caus], gt[caus])[0, 1]
        print(f"  {name:9s}: causal corr with truth = {r:.3f}, "
              f"acausal max|beta| = {np.max(np.abs(k[~caus])):.3f}")

    # 7. figure ---------------------------------------------------------------
    t_rel = -lag_s
    o = np.argsort(t_rel); tx = t_rel[o]
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    ax.axhline(0, color="k", lw=0.8)
    ax.plot(tx, scale_to(kernel_p)[o], "--", color="#d9331e", lw=2.2,
            label=(r"ground truth $e^{-\tau/%g}$" % TAU_GEN) if not np.isinf(TAU_GEN)
                  else "ground truth (flat)")
    ax.plot(tx, kernel_p[o], "-", color="0.1", lw=2.6, label=f"P-spline curvature (λ={lam_p:.0f})")
    ax.plot(tx, kernel_l2[o], "-", color="#1c6bae", lw=2.0, label=f"sklearn ridge L2 (C={C_l2:.2g})")
    ax.plot(tx, kernel_l1[o], "-", color="#2aa36b", lw=2.0, label=f"sklearn lasso L1 ({nnz}/{K} nz)")
    ax.set_xlabel("time relative to un-freeze (s)")
    ax.set_ylabel(r"$\beta(\tau)$")
    ax.set_title("recovered kernel: sklearn built-in penalties vs curvature penalty")
    ax.legend(frameon=False, loc="best")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = "kernel_recovery_sklearn_py.png"
    fig.savefig(out, dpi=150)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
