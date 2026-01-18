function out = ddm_pdf_from_trace(params, scaled_drift, fixed)
% params = [lambda, bound, ndt]
% scaled_drift: the stored trace (column vector, coarse dt)
% fixed: struct containing dt, T_max, T_trunc, etc.

lambda      = params(1);
bound       = params(2);
ndt         = params(3); % Non-decision time in seconds

% --- 1. PDE Grid Setup ---
grid_size = round((bound - fixed.x_min) / fixed.dx) + 1;
start_idx = round((fixed.x0    - fixed.x_min) / fixed.dx) + 1;

mex_p = [fixed.dt, fixed.dx, fixed.sigma_sq, fixed.x0, ...
         fixed.x_min, grid_size, start_idx, lambda];

% --- 2. Solve PDE (decision time only) ---
drift_used = scaled_drift(:);
[p_dist, survival_at_Tmax] = leaky_pde_robust(drift_used, mex_p);

% Decision-time grid
t_dec = (0:fixed.dt:fixed.T_max)';

% --- 3. Apply Non-Decision Time (NDT) ---
ndt_idx = round(ndt / fixed.dt);

if ndt_idx < 0
    error('Negative NDT is not physically meaningful.');
end

% Shift PDF safely (no silent mass loss)
pdf_shifted = zeros(size(p_dist));
if ndt_idx < length(p_dist)
    pdf_shifted(ndt_idx+1:end) = p_dist(1:end-ndt_idx);
end

% --- 4. Survival after NDT (critical fix) ---
% Observed survival at T_max:
% P(T_obs > T_max) = P(T_dec > T_max - ndt)
if ndt_idx < length(p_dist)
    overflow_mass = sum(p_dist(end-(ndt_idx-1):end)) * fixed.dt;
    survival_ndt  = survival_at_Tmax + overflow_mass;
else
    survival_ndt = 1.0;
end

% --- 5. Truncation Renormalization (in observed time) ---
idx_trunc_obs = round(fixed.T_trunc / fixed.dt);

if idx_trunc_obs > 1
    prob_missing = sum(pdf_shifted(1:idx_trunc_obs)) * fixed.dt;
else
    prob_missing = 0;
end

norm_factor = max(1e-12, 1.0 - prob_missing);

pdf_renorm  = pdf_shifted / norm_factor;
surv_renorm = survival_ndt / norm_factor;

% --- 6. Package Output ---
out.t            = t_dec;          % observed time grid
out.pdf_raw      = p_dist(:);            % decision-time PDF
out.pdf_shifted  = pdf_shifted;          % NDT-shifted (unnormalized)
out.pdf          = pdf_renorm;            % final likelihood PDF
out.survival_raw = survival_at_Tmax;     % decision-time survival
out.survival     = surv_renorm;           % observed-time survival
out.norm_factor  = norm_factor;
out.ndt_idx      = ndt_idx;

% --- 7. Discrete PMF (for likelihood comparison) ---
out.pmf = out.pdf * fixed.dt;
out.pmf = max(out.pmf, 1e-12);

end
