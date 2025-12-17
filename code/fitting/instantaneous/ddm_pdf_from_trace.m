function out = ddm_pdf_from_trace(params, drift_coarse, fixed)
% params = [drift_scale, lambda, bound]
% drift_coarse: the stored trace in drifts_cell{i} (column vector, coarse dt)
% fixed: your struct

drift_scale = params(1);
lambda      = params(2);
bound       = params(3);

% --- Dynamic grid for THIS bound ---
grid_size = round((bound - fixed.x_min) / fixed.dx) + 1;
start_idx = round((fixed.x0   - fixed.x_min) / fixed.dx) + 1;

mex_p = [fixed.dt, fixed.dx, fixed.sigma_sq, fixed.x0, ...
         fixed.x_min, grid_size, start_idx, lambda];

% --- Drift seen by PDE ---
drift_used = drift_coarse(:) * drift_scale;

% --- Solve PDE ---
[p_dist, survival_at_Tmax] = leaky_pde_robust(drift_used, mex_p);

t = (0:fixed.dt:fixed.T_max)';

% --- Apply the SAME truncation renormalization as in your likelihood ---
idx_trunc = round(fixed.T_trunc / fixed.dt);
if idx_trunc > 1
    prob_missing = sum(p_dist(1:idx_trunc)) * fixed.dt;
else
    prob_missing = 0;
end
norm_factor = max(1e-12, 1.0 - prob_missing);

pdf_renorm = p_dist(:) / norm_factor;
surv_renorm = survival_at_Tmax / norm_factor;

% Package
out.t = t;
out.pdf_raw = p_dist(:);
out.pdf = pdf_renorm;
out.survival_raw = survival_at_Tmax;
out.survival = surv_renorm;
out.norm_factor = norm_factor;
end
