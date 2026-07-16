"""
recover_demo.py — ground-truth recovery of the social-motion hazard kernel.

Self-contained NumPy/Matplotlib mirror of kernel_recovery_ddm.m:

  1. Make autocorrelated social-motion traces (stand-ins for the real per-fly sm).
  2. Simulate un-freeze times from a leaky accumulator (DDM) with FIXED params,
     driven by that motion -> a KNOWN kernel beta_gt(tau) ~ exp(-tau/tau_gen).
  3. Recover beta(tau) with the person-period logistic GLM (raised-cosine basis,
     2nd-difference smoothness penalty, bout-grouped CV for lambda).
  4. Overlay recovered vs ground-truth kernel and the CV-vs-lambda curve.

Run:  python recover_demo.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import hazard_kernel as hk


# ── config ───────────────────────────────────────────────────────────────────
FPS = 60
KERNEL_PAST_S = 4.0
KERNEL_FUTURE_S = 1.0        # acausal control -> should recover ~0
DT_FRAMES = 6
ENTRY_FR = 30                # left-truncation (0.5 s)
NB_BASIS = 12
BUMP_WIDTH = 1.25
NB_BASELINE = 8
CV_FOLDS = 5
GRID_ANCHOR = "offset"
MASK_PREONSET = True

# generative DDM (fixed, known) — tau_gen = np.inf -> perfect DDM (flat kernel)
TAU_GEN = 1.0
MU0, BETA, SIGMA, THETA = 0.0, 15.0, 0.4, 2.0   # low noise + strong social drive
DT = 1.0 / FPS
MAXT_FR = int(round(10.5 * FPS))

N_BOUTS = 2000               # distinct real-fly templates
N_AUG = 5                    # sim trajectories per template (more events)
SEED = 7


def make_social_motion(n_bouts, rng):
    """Autocorrelated, non-negative per-bout motion traces (onset forward)."""
    sm_list = []
    for _ in range(n_bouts):
        length = int(rng.integers(180, MAXT_FR))
        # smoothed AR(1)-ish positive signal
        z = rng.standard_normal(length)
        a = 0.94
        x = np.zeros(length)
        for k in range(1, length):
            x[k] = a * x[k - 1] + np.sqrt(1 - a * a) * z[k]
        sm = np.abs(x) + 0.3 * rng.random()          # keep >= 0, small per-bout offset
        sm_list.append(sm.astype(float))
    return sm_list


def main():
    rng = np.random.default_rng(SEED)

    # 1. social motion (templates) --------------------------------------------
    sm_list_templates = make_social_motion(N_BOUTS, rng)

    # augment: replicate each template N_AUG times (shared bout id for grouped CV)
    sm_list = []
    bout_of_row = []
    for rep in range(N_AUG):
        for b in range(N_BOUTS):
            sm_list.append(sm_list_templates[b])
            bout_of_row.append(b)          # shared template id across reps
    n_sim = len(sm_list)

    # 2. simulate DDM un-freeze times -----------------------------------------
    lam = 0.0 if np.isinf(TAU_GEN) else 1.0 / TAU_GEN
    drives = [MU0 + BETA * (sm / 10.0) for sm in sm_list]     # per-frame drift rate
    durations, censored = hk.simulate_accumulator(
        drives, THETA, SIGMA, lam, DT, seed=SEED + 100)
    print(f"Simulated {n_sim} bouts | censored {100*censored.mean():.1f}% | "
          f"median dur {np.median(durations[~censored])/FPS:.2f}s")

    # 3. lag grid + global z-scoring of sm ------------------------------------
    back = int(round(KERNEL_PAST_S * FPS))
    fwd = int(round(KERNEL_FUTURE_S * FPS))
    lag_fr = np.arange(-fwd, back + 1)          # + = past, - = future
    lag_s = lag_fr / FPS
    caus = lag_fr >= 0
    allsm = np.concatenate(sm_list_templates)
    sm_mu, sm_sd = allsm.mean(), allsm.std()

    # 4. build the person-period design ---------------------------------------
    S, y, tinf, gid_local = hk.build_design(
        durations, censored, sm_list, lag_fr, FPS, DT_FRAMES, ENTRY_FR,
        GRID_ANCHOR, MASK_PREONSET, sm_mu, sm_sd)
    gid = np.array([bout_of_row[i] for i in gid_local])       # template id (grouped CV)
    print(f"Design: {S.shape[0]} at-risk intervals, {int(y.sum())} events "
          f"(hazard {y.mean():.3f})")

    # 5. basis, projection, baseline nuisance ---------------------------------
    B = hk.raised_cosine_basis(lag_fr, NB_BASIS, BUMP_WIDTH)
    K = B.shape[1]
    Xk = S @ B
    # baseline hazard over time-in-freeze (rank-normalised), drop first col
    ut = (np.argsort(np.argsort(tinf)) + 1) / tinf.size
    Bb = hk.raised_cosine_basis(ut * (back + fwd) - fwd, NB_BASELINE, 1.25)[:, 1:]
    X = np.column_stack([np.ones(S.shape[0]), Xk, Bb])
    ker = np.arange(1, 1 + K)                                  # kernel coefficient cols

    # penalty on the kernel block only
    D = hk.second_difference(K)
    Pk = np.zeros((X.shape[1], X.shape[1]))
    Pk[np.ix_(ker, ker)] = D.T @ D

    # 6. choose lambda by grouped CV, then fit --------------------------------
    lam_grid = np.logspace(-1, 5, 10)
    lam, cvll = hk.select_lambda(X, y, Pk, gid, lam_grid, CV_FOLDS, seed=SEED, rule="argmax")
    b, V = hk.penalized_logit(X, y, lam * Pk)
    bk = b[ker]
    kernel = B @ bk
    kse = np.sqrt(np.sum((B @ V[np.ix_(ker, ker)]) * B, axis=1))
    print(f"Selected lambda = {lam:.3g}")

    # 7. ground truth + recovered timescale -----------------------------------
    gt = np.zeros_like(kernel)
    gt[caus] = 1.0 if np.isinf(TAU_GEN) else np.exp(-lag_s[caus] / TAU_GEN)
    a_scale = (gt[caus] @ kernel[caus]) / max(gt[caus] @ gt[caus], 1e-12)
    gt_scaled = a_scale * gt
    tau_hat = hk.fit_leak_tau(kernel, lag_s, caus)
    print(f"tau_gen = {TAU_GEN}  |  tau_hat = {tau_hat:.2f}s  |  "
          f"peak = {np.max(np.abs(kernel[caus])):.3f}  |  "
          f"acausal max|beta| = {np.max(np.abs(kernel[~caus])):.3f}")

    # 8. figure ----------------------------------------------------------------
    t_rel = -lag_s
    order = np.argsort(t_rel)
    tx = t_rel[order]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 4.6))

    axL.semilogx(lam_grid, cvll, "-o", color="#1c6bae", lw=2)
    axL.axvline(lam, ls=":", color="#d95f1e", lw=1.5)
    axL.set_xlabel("smoothness  lambda")
    axL.set_ylabel("held-out log-likelihood / obs")
    axL.set_title("bout-grouped CV")
    axL.spines[["top", "right"]].set_visible(False)

    axR.fill_between(tx, (kernel + 1.96 * kse)[order], (kernel - 1.96 * kse)[order],
                     color="0.2", alpha=0.15, lw=0)
    axR.axhline(0, color="k", lw=0.8)
    axR.plot(tx, gt_scaled[order], "--", color="#d9331e", lw=2.2,
             label=("ground truth (flat)" if np.isinf(TAU_GEN)
                    else fr"ground truth $e^{{-\tau/{TAU_GEN:g}}}$"))
    axR.plot(tx, kernel[order], "-", color="0.1", lw=2.6,
             label=fr"recovered ($\hat\tau$={tau_hat:.2f}s)")
    axR.set_xlabel("time relative to un-freeze (s)")
    axR.set_ylabel(r"$\beta(\tau)$")
    axR.set_title("recovered vs ground-truth kernel")
    axR.legend(frameon=False, loc="best")
    axR.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    out = "kernel_recovery_ddm_py.png"
    fig.savefig(out, dpi=150)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
