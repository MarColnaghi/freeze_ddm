function lik = ddm_freeze_loglik_from_trace_ultrafast( ...
        params, scaled_drift, fixed, rt)
% params = [lambda, bound, ndt]
% rt     = observed RT (NaN or >= T_max => censored)

lambda = params(1);
bound  = params(2);
ndt    = params(3);

dt   = fixed.dt;
Tmax = fixed.T_max;

%% --- 1. PDE grid ---
grid_size = round((bound - fixed.x_min) / fixed.dx) + 1;
start_idx = round((fixed.x0 - fixed.x_min) / fixed.dx) + 1;

mex_p = [dt, fixed.dx, fixed.sigma_sq, fixed.x0, ...
         fixed.x_min, grid_size, start_idx, lambda];

%% --- 2. Solve PDE (decision time only) ---
p_dist = leaky_pde_robust(scaled_drift(:), mex_p);
Nt = length(p_dist);

% Precompute CDF ONCE
cdf = cumsum(p_dist) * dt;

%% --- 3. Truncation normalization (VBMC-consistent) ---
trunc_dec_time = fixed.T_trunc - ndt;

if trunc_dec_time <= 0
    missing_mass = 0;
else
    idx_trunc = floor(trunc_dec_time / dt);
    if idx_trunc < 1
        missing_mass = 0;
    elseif idx_trunc >= Nt
        missing_mass = 1.0;
    else
        missing_mass = cdf(idx_trunc);
    end
end

Z = max(1e-12, 1.0 - missing_mass);

%% --- 4. Likelihood ---
% ---- Censored ----
if isnan(rt) || rt >= (Tmax - dt)

    cutoff_dec = Tmax - ndt;

    if cutoff_dec <= 0
        surv = 1.0;
    else
        idx_cut = floor(cutoff_dec / dt);
        if idx_cut < 1
            absorbed = 0.0;
        elseif idx_cut >= Nt
            absorbed = 1.0;
        else
            absorbed = cdf(idx_cut);
        end
        surv = 1.0 - absorbed;
    end

    lik = surv / Z;

% ---- Observed RT ----
else
    decision_time = rt - ndt;

    if decision_time <= 0
        lik = 0;
    else
        idx = floor(decision_time / dt);
        if idx < 1 || idx >= Nt
            lik = 0;
        else
            lik = p_dist(idx) / Z;
        end
    end
end

end
