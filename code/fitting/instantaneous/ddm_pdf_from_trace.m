function out = ddm_pdf_from_trace(params, scaled_drift, fixed)
% params = [lambda, bound, ndt]
% scaled_drift : drift trace (column vector, coarse dt)
% fixed        : struct with fields:
%                dt, dx, sigma_sq, x0, x_min,
%                T_max, T_trunc

lambda = params(1);
bound  = params(2);
ndt    = params(3);   % non-decision time (seconds)

dt = fixed.dt;

%% --- 1. PDE Grid Setup ---
grid_size = round((bound - fixed.x_min) / fixed.dx) + 1;
start_idx = round((fixed.x0   - fixed.x_min) / fixed.dx) + 1;

mex_p = [fixed.dt, fixed.dx, fixed.sigma_sq, fixed.x0, ...
         fixed.x_min, grid_size, start_idx, lambda];

%% --- 2. Solve PDE (decision time only) ---
drift_used = scaled_drift(:);

[p_dist, survival_at_Tmax] = leaky_pde_robust(drift_used, mex_p);

% Decision-time grid
t_dec = (0:dt:fixed.T_max)';

%% --- 3. Apply Non-Decision Time (NDT) ---
ndt_idx = round(ndt / dt);
if ndt_idx < 0
    error('Negative NDT is not physically meaningful.');
end

pdf_shifted = zeros(size(p_dist));
if ndt_idx < length(p_dist)
    pdf_shifted(ndt_idx+1:end) = p_dist(1:end-ndt_idx);
end

%% --- 4. Survival after NDT ---
% P(T_obs > T_max) = P(T_dec > T_max - ndt)

if ndt_idx < length(p_dist)
    overflow_mass = sum(p_dist(end-(ndt_idx-1):end)) * dt;
    survival_ndt  = survival_at_Tmax + overflow_mass;
else
    survival_ndt = 1.0;
end

%% --- 5. Truncation conditioning (CORRECT) ---
% Condition on T_obs > T_trunc

idx_trunc = floor(fixed.T_trunc / dt) + 1;

if idx_trunc <= length(pdf_shifted)
    mass_after_trunc = sum(pdf_shifted(idx_trunc:end)) * dt;
else
    mass_after_trunc = 0;
end

survival_trunc = mass_after_trunc + survival_ndt;
norm_factor    = max(1e-12, survival_trunc);

pdf_renorm  = pdf_shifted / norm_factor;
surv_renorm = survival_ndt  / norm_factor;

%% --- 6. Package Output ---
out.t            = t_dec;                % observed-time grid
out.pdf_raw      = p_dist(:);            % decision-time PDF
out.pdf_shifted  = pdf_shifted;          % NDT-shifted (unnormalized)
out.pdf          = pdf_renorm;            % truncated likelihood PDF
out.survival_raw = survival_at_Tmax;     % decision-time survival
out.survival     = surv_renorm;           % observed-time survival
out.norm_factor  = norm_factor;
out.ndt_idx      = ndt_idx;
out.idx_trunc    = idx_trunc;

%% --- 7. Discrete PMF (for likelihoods) ---
out.pmf = out.pdf * dt;
out.pmf = max(out.pmf, 1e-12);

end