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
% NB: grid is recomputed from the (possibly free / fitted) bound here, which
% intentionally OVERRIDES fixed.grid_size / fixed.start_idx from param_res.
grid_size = round((bound - fixed.x_min) / fixed.dx) + 1;
start_idx = round((fixed.x0   - fixed.x_min) / fixed.dx) + 1;

mex_p = [fixed.dt, fixed.dx, fixed.sigma_sq, fixed.x0, ...
         fixed.x_min, grid_size, start_idx, lambda];

%% --- 2. Solve PDE (decision time only) ---
drift_used = scaled_drift(:);

[p_dist, survival_at_Tmax] = leaky_pde_robust(drift_used, mex_p);

% Decision-time grid (aligned 1:1 with the PDF / drift trace, so out.t and
% out.pdf ALWAYS have equal length -- callers must not slice one to fit).
t_dec = (0:numel(p_dist)-1)' * dt;

%% --- 3. Apply Non-Decision Time (NDT) ---
% Sub-bin (fractional) shift: linearly blend the floor/ceil whole-bin shifts so
% the likelihood is SMOOTH in ndt (helps gradient/simplex fitting). For an
% integer ndt/dt this reduces exactly to the old round-to-nearest shift.
if ndt < 0
    error('Negative NDT is not physically meaningful.');
end
ndt_bins = ndt / dt;
n0 = floor(ndt_bins);
w  = ndt_bins - n0;            % fractional part in [0,1)
np = length(p_dist);

pdf_shifted = zeros(size(p_dist));
if n0 < np
    pdf_shifted(n0+1:end) = pdf_shifted(n0+1:end) + (1 - w) * p_dist(1:np-n0);
end
if (n0 + 1) < np
    pdf_shifted(n0+2:end) = pdf_shifted(n0+2:end) + w * p_dist(1:np-n0-1);
end

%% --- 4. Survival after NDT ---
% Decision mass shifted PAST T_max becomes censored. This is just the mass that
% no longer fits on the [0, T_max] grid after the shift -- valid for fractional
% shifts too, and reduces to summing the last ndt bins for an integer shift.
overflow_mass = (sum(p_dist) - sum(pdf_shifted)) * dt;
survival_ndt  = survival_at_Tmax + max(0, overflow_mass);

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

% The conditional density is defined only for T_obs > T_trunc; zero the earlier
% bins so out.pdf is a PROPER truncated density (its integral over the window
% plus surv_renorm = 1). Callers no longer need to slice at idx_trunc.
if idx_trunc > 1
    pdf_renorm(1:idx_trunc-1) = 0;
end

%% --- 6. Package Output ---
out.t            = t_dec;                % observed-time grid
out.pdf_raw      = p_dist(:);            % decision-time PDF
out.pdf_shifted  = pdf_shifted;          % NDT-shifted (unnormalized)
out.pdf          = pdf_renorm;            % truncated likelihood PDF
out.survival_raw = survival_at_Tmax;     % decision-time survival
out.survival     = surv_renorm;           % observed-time survival
out.norm_factor  = norm_factor;
out.ndt_idx      = n0;                    % integer part of the ndt shift (bins)
out.ndt_frac     = w;                     % sub-bin fractional part
out.idx_trunc    = idx_trunc;

%% --- 7. Discrete PMF (for likelihoods) ---
out.pmf = out.pdf * dt;
out.pmf = max(out.pmf, 1e-12);            % floor for log-likelihood safety
if idx_trunc > 1
    out.pmf(1:idx_trunc-1) = 0;           % below-truncation region is impossible
end

end