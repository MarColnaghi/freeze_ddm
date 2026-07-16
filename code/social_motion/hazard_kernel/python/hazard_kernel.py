"""
hazard_kernel.py — a compact NumPy port of the social-motion hazard-kernel GLM.

Discrete-time survival / person-period logistic model for the un-freeze hazard,
driven by a temporal kernel over recent social motion:

    logit P(break at interval i | still frozen) = b0 + sum_tau beta(tau) * sm(t - tau)

The kernel beta(tau) is expanded in a small raised-cosine basis (beta = B @ b),
fit by penalised IRLS (Newton) with an optional 2nd-difference smoothness
penalty lambda * D'D, and lambda is chosen by bout-grouped cross-validation.

Only numpy is required here (matplotlib is used by the demo script).
"""

import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# Basis + penalty
# ─────────────────────────────────────────────────────────────────────────────
def raised_cosine_basis(lag_fr, n_basis, width_factor=1.25):
    """Evenly spaced, overlapping raised-cosine (Hann) bumps over the lag axis.

    Returns B of shape (nLag, n_basis); each column is one bump, peak 1.
    Adjacent bumps overlap (width_factor > 1), which is what makes beta = B@b smooth.
    """
    lag_fr = np.asarray(lag_fr, float)
    centres = np.linspace(lag_fr.min(), lag_fr.max(), n_basis)
    spacing = centres[1] - centres[0] if n_basis > 1 else (lag_fr.ptp() + 1)
    width = width_factor * spacing
    B = np.zeros((lag_fr.size, n_basis))
    for j, c in enumerate(centres):
        d = np.abs(lag_fr - c)
        m = d < width
        B[m, j] = 0.5 * (1 + np.cos(np.pi * (lag_fr[m] - c) / width))
    return B


def second_difference(K):
    """(K-2) x K second-difference operator D (rows [1 -2 1]); ||D b||^2 = roughness."""
    return np.diff(np.eye(K), n=2, axis=0)


# ─────────────────────────────────────────────────────────────────────────────
# Penalised logistic regression by IRLS (Newton), with covariance
# ─────────────────────────────────────────────────────────────────────────────
def penalized_logit(X, y, P, maxit=100, tol=1e-8, ridge=1e-8):
    """MAP logistic regression: maximise loglik(b) - 0.5 * b' P b.

    Each iteration is one weighted least-squares solve (IRLS):
        w = h(1-h);  H = X'WX + P;  g = X'(y-h) - P b;  b <- b + H^{-1} g
    Returns (b, V) where V = inv(H) at the solution (coefficient covariance).
    """
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    n, p = X.shape
    b = np.zeros(p)
    R = ridge * np.eye(p)
    H = None
    for _ in range(maxit):
        eta = np.clip(X @ b, -30, 30)
        mu = 1.0 / (1.0 + np.exp(-eta))
        w = np.maximum(mu * (1 - mu), 1e-9)
        H = (X * w[:, None]).T @ X + P + R
        g = X.T @ (y - mu) - P @ b
        step = np.linalg.solve(H, g)
        b = b + step
        if np.max(np.abs(step)) < tol:
            break
    V = np.linalg.inv(H)
    return b, V


def _bernoulli_loglik_per_obs(X, y, b):
    eta = np.clip(X @ b, -30, 30)
    mu = np.clip(1.0 / (1.0 + np.exp(-eta)), 1e-9, 1 - 1e-9)
    return np.mean(y * np.log(mu) + (1 - y) * np.log(1 - mu))


# ─────────────────────────────────────────────────────────────────────────────
# Bout-grouped cross-validation to pick lambda
# ─────────────────────────────────────────────────────────────────────────────
def grouped_cv_loglik(X, y, P, groups, k_folds, seed=0):
    """Held-out Bernoulli log-likelihood/obs, folds grouped by `groups` (bout id)."""
    rng = np.random.default_rng(seed)
    uniq = np.unique(groups)
    fold_of_group = rng.integers(0, k_folds, size=uniq.size)
    g2f = dict(zip(uniq.tolist(), fold_of_group.tolist()))
    rowfold = np.array([g2f[g] for g in groups])
    lls = []
    for k in range(k_folds):
        te = rowfold == k
        tr = ~te
        if not te.any() or not tr.any():
            continue
        b, _ = penalized_logit(X[tr], y[tr], P)
        lls.append(_bernoulli_loglik_per_obs(X[te], y[te], b))
    return float(np.mean(lls)) if lls else np.nan


def select_lambda(X, y, Pk, groups, lam_grid, k_folds=5, seed=0, rule="argmax"):
    """Return (lambda, cvll_array). `Pk` is the *unit* penalty; scaled by each lambda.

    rule='argmax' picks the best CV score; rule='1se' picks the largest lambda
    within one standard error of the best (sparser / smoother, robust when flat).
    """
    cvll = np.array([grouped_cv_loglik(X, y, lam * Pk, groups, k_folds, seed)
                     for lam in lam_grid])
    ibest = int(np.nanargmax(cvll))
    if rule == "1se":
        # crude 1-SE using spread across the grid neighbourhood as a proxy SE
        se = np.nanstd(cvll) / np.sqrt(max(k_folds, 1))
        ok = np.where(cvll >= cvll[ibest] - se)[0]
        isel = ok.max()  # largest lambda within 1 SE
    else:
        isel = ibest
    return lam_grid[isel], cvll


# ─────────────────────────────────────────────────────────────────────────────
# Leaky-accumulator (DDM) ground-truth simulator — vectorised across bouts
# ─────────────────────────────────────────────────────────────────────────────
def simulate_accumulator(drives, theta, sigma, lam, dt, seed=0):
    """First-passage times of a leaky accumulator integrating each bout's drive.

        x(k+1) = x(k)(1 - lam*dt) + drive(k)*dt + sigma*sqrt(dt)*randn
        un-freeze = first k with x >= theta   (else censored at that bout's length)

    lam = 0     -> perfect DDM (flat kernel);
    lam = 1/tau -> leaky (exponential kernel with time-constant tau).

    `drives` is a list of 1-D arrays (per-bout, per-frame drift RATE).
    Returns (durations_1based, censored_bool).
    """
    rng = np.random.default_rng(seed)
    n = len(drives)
    lengths = np.array([len(d) for d in drives])
    T = int(lengths.max())
    D = np.zeros((n, T))
    for b, d in enumerate(drives):
        D[b, : len(d)] = d
    decay = 1 - lam * dt
    nscale = sigma * np.sqrt(dt)
    x = np.zeros(n)
    rt = np.full(n, -1, dtype=int)
    for k in range(T):
        active = (rt < 0) & (k < lengths)
        noise = rng.standard_normal(n)
        xnew = x * decay + D[:, k] * dt + nscale * noise
        x = np.where(active, xnew, x)
        newly = active & (x >= theta)
        rt[newly] = k + 1
    durations = np.where(rt > 0, rt, lengths)
    censored = rt < 0
    return durations, censored


# ─────────────────────────────────────────────────────────────────────────────
# Person-period design from (synthetic) durations + real/sim social motion
# ─────────────────────────────────────────────────────────────────────────────
def build_design(durations, censored, sm_list, lag_fr, fps, dt_frames,
                 entry_fr, grid_anchor="offset", mask_preonset=True,
                 sm_mu=0.0, sm_sd=1.0):
    """Stack every at-risk interval of every bout into the long design.

    Each bout's `sm_list[b]` is the social motion from freeze onset forward
    (index 0 = onset = frame 1). `lag_fr` is positive for the past, negative for
    the (acausal) future. Returns S (rows x nLag), y (0/1), tinf (s), gid (bout).
    """
    lag_fr = np.asarray(lag_fr, int)
    nLag = lag_fr.size
    rows_S, y, tinf, gid = [], [], [], []
    for b in range(len(sm_list)):
        sm = np.asarray(sm_list[b], float)
        L = sm.size
        dur = int(durations[b])
        if dur < entry_fr:
            continue
        base = entry_fr if grid_anchor == "entry" else dur
        kk = np.arange(np.ceil((1 - base) / dt_frames),
                       np.floor((dur - base) / dt_frames) + 1)
        allg = base + kk * dt_frames
        grid = allg[allg >= entry_fr]
        if grid_anchor == "entry":
            if grid.size == 0 or grid[-1] < dur:
                grid = np.append(grid, dur)
        else:
            if grid.size == 0 or grid[0] > entry_fr:
                grid = np.insert(grid, 0, entry_fr)
        grid = grid.astype(int)
        for gi, f in enumerate(grid):
            idx = f - lag_fr  # absolute 1-based frame per lag (onset = frame 1)
            if idx.max() > L:              # missing FUTURE data -> can't build this interval
                continue
            pre = idx < 1                  # lags reaching before onset (no data available)
            if pre.any() and not mask_preonset:
                continue                   # would need real pre-onset data this port doesn't have
            s = np.zeros(nLag)
            ok = ~pre
            s[ok] = (sm[idx[ok] - 1] - sm_mu) / sm_sd
            if mask_preonset:
                s[lag_fr > f - 1] = 0.0    # (== pre) neutral 0 = global mean
            if np.isnan(s).any():
                continue
            rows_S.append(s)
            tinf.append(f / fps)
            gid.append(b)
            y.append(1 if (gi == len(grid) - 1 and not censored[b]) else 0)
    return (np.asarray(rows_S), np.asarray(y, float),
            np.asarray(tinf), np.asarray(gid, int))


# ─────────────────────────────────────────────────────────────────────────────
# Utilities
# ─────────────────────────────────────────────────────────────────────────────
def fit_leak_tau(kernel, lag_s, caus):
    """Coarse exponential fit beta(tau) ~ A*exp(-tau/tau)+c on causal lags -> tau."""
    t = lag_s[caus]
    kc = kernel[caus]
    best_tau, best_sse = np.nan, np.inf
    for tau in np.logspace(np.log10(0.08), np.log10(8.0), 80):
        Xf = np.column_stack([np.exp(-t / tau), np.ones_like(t)])
        coef, *_ = np.linalg.lstsq(Xf, kc, rcond=None)
        sse = np.sum((kc - Xf @ coef) ** 2)
        if sse < best_sse:
            best_sse, best_tau = sse, tau
    return best_tau
