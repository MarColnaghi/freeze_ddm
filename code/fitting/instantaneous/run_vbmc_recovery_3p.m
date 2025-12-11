% Check for VBMC
if isempty(which('vbmc'))
    error('VBMC is not in the path. Please install/add it.');
end

%% 1. Setup Parameters
rng(randi([1 60]));

% --- GROUND TRUTH ---
% We will fit:
% 1. Mean Drift (base level before noise)
% 2. Lambda (Leak)
true_drift_scale = 1.1;
true_lambda = .20;
true_bound = 0.6;


% -- Likelihood / Signal Parameters (Coarse) --
fixed.dt = 0.004;       % 1ms for signal and PDE
fixed.dx = 0.01;
fixed.sigma_sq = 1.0;
%     fixed.bound = 1.0;
fixed.x0 = 0.0;

% -- Simulation Parameters (Fine) --
sim.dt = 0.00005;       % 0.05ms for ground truth simulation

% PDE Grid Setup
x_min = -7;
grid_size = round((true_bound - x_min) / fixed.dx) + 1;
start_idx = round((fixed.x0 - x_min) / fixed.dx) + 1;

% Save for objective function
fixed.x_min = x_min;
fixed.grid_size = grid_size;
fixed.start_idx = start_idx;

% Truncation and T_max
fixed.T_trunc = 0;
fixed.T_max = 10.0;


fprintf('--- Ground Truth ---\nMean Drift: %.2f\nLambda: %.2f\n', true_drift_scale, true_lambda);

%% 2. Generate Synthetic Data (Parallel Data Gen)
N_trials = 5000;

% Time vectors
t_coarse = 0:fixed.dt:fixed.T_max;
t_fine   = 0:sim.dt:fixed.T_max;

fprintf('\nGenerating %d trials with complex drift...\n', N_trials);
fprintf('  - Signal dt: %.5fs\n  - Sim dt:    %.5fs\n', fixed.dt, sim.dt);

drifts_cell = cell(N_trials, 1); % Stores the COARSE signal (what the model sees)
true_rt = zeros(N_trials, 1);

% Params vector for Sim: [dt_fine, sigma, x0, theta, lambda]
sim_params_vec = [sim.dt, fixed.sigma_sq, fixed.x0, true_bound, true_lambda];

tic;
%     base_signal_zero = generate_complex_drift(t_coarse, fixed.dt, 0.0);
for i = 1:N_trials
    % 1. Generate the COARSE signal using your specific logic
    % Generate signal with ZERO mean
    base_signal_zero = generate_complex_drift(t_coarse, fixed.dt, 1.0);

    % Store this pure noise trace
    drifts_cell{i} = base_signal_zero;

    % Add TRUE mean for the Ground Truth Simulation
    sim_signal_coarse = base_signal_zero * true_drift_scale;

    % Interpolate THIS to fine
    fine_signal = interp1(t_coarse(:), sim_signal_coarse(:), t_fine(:), 'nearest');

    % Handle NaNs at the end if dimensions don't match perfectly
    fine_signal(isnan(fine_signal)) = sim_signal_coarse(end);

    % 3. Run High-Res Simulation (C++ MEX)
    % We use the FINE signal and FINE dt
    % CRITICAL FIX: Pass a unique seed for this trial!
    % Using i (loop index) ensures it's unique but reproducible.
    % Adding a large offset prevents overlap with other RNGs if used elsewhere.
    unique_seed = uint64(randi([1 60]) + i * 1000);

    res = sim_ddm_seeded(fine_signal, sim_params_vec, 1, unique_seed);

    if isnan(res)
        true_rt(i) = fixed.T_max;
    else
        true_rt(i) = res;
    end
end
fprintf('Data generation complete in %.2f seconds.\n', toc);

%% 3. Setup VBMC

% Order: [Drift, Lambda, Bound]
% Note: Bound should be positive.

%       Drift   Lambda   Bound
LB  = [ 0.01,   0,     0.3 ];
UB  = [ 5.0,    5.0,     4.0 ];

% Plausible
PLB = [ 0.25,   +.05,    0.5 ];
PUB = [ 4.0,    4.0,     2.5 ];

x0_start = [1.5, .5, 1.0];

% Update Handles
%     log_prior_fun = @(p) log_correlated_prior_3D_robust(p, LB, UB, PLB, PUB);
log_prior_fun = @(p) log_trapezoidal_prior(p, LB, UB, PLB, PUB);

% Note: Likelihood function now needs to handle the dynamic bound internally
log_lik_fun = @(p) log_joint_likelihood_3D(p, drifts_cell, true_rt, fixed);

% C. Define The Posterior (Target for VBMC)
% Posterior = Likelihood + Prior
log_posterior_fun = @(p) log_lik_fun(p) + log_prior_fun(p);

options = vbmc('defaults');
options.Display = 'iter';
options.Plot = false;
options.SpecifyTargetNoise = false;
options.TolStableCount = 80;
options.MinFinalComponents = 50;

% 1. Force a high-resolution initial map
options.FunEvalStart = 50;
% 2. Force it to run longer (don't stop on stability early)
options.MinFunEvals = 250;
options.MaxFunEvals = 800;
% 3. Force higher precision before declaring convergence
options.TolImprovement = 0.002;
% 4. Start with a more complex mixture model (smoother tails)
options.Kwarmup = 5;


%% 4. Run Inference
fprintf('\nStarting VBMC Inference...\n');
tic;
[vp, elbo, ~] = vbmc(log_posterior_fun, x0_start, LB, UB, PLB, PUB, options);
total_time = toc;

fprintf('\nInference complete in %.2f minutes.\n', total_time/60);

%% 5. Results
[post_mean,post_cov] = vbmc_moments(vp);
post_std = sqrt(diag(post_cov));

fprintf('\n--- Results ---\n');
fprintf('True Drift Scale: %.2f | Est: %.2f (+/- %.2f)\n', true_drift_scale, post_mean(1), post_std(1));
fprintf('True Lambda:      %.2f | Est: %.2f (+/- %.2f)\n', true_lambda, post_mean(2), post_std(2));
fprintf('True Bound:       %.2f | Est: %.2f (+/- %.2f)\n', true_bound, post_mean(3), post_std(3));

vbmc_plot(vp); sgtitle('Parameter Recovery');

%obtained log likelihood
obtained_LL = log_joint_likelihood_3D(post_mean, drifts_cell, true_rt, fixed);
%true log likelihood
true_LL = log_joint_likelihood_3D([true_drift_scale,true_lambda,true_bound], drifts_cell, true_rt, fixed);
fprintf('Obtained Log Likelihood: %.2f\n', obtained_LL);
fprintf('True Log Likelihood:     %.2f\n', true_LL);

%% plot drift

figure('color','w'); hold on
plot(t_coarse, drifts_cell{1}, 'k-', 'LineWidth', 1);
title('Time-Varying Drift');
xlabel('Time (s)');
ylabel('Drift (\mu)');

res = 20;
drifts = linspace(1.25, 1.35, res);
bounds = linspace(2.07, 2.13, res);
LL = nan(res,res);
i = 0;
for idx_drift_scale = 1:res
    for idx_bound = 1:res
        i = i + 1;
        fprintf('done %d \n', i)
        drift_scale = drifts(idx_drift_scale);
        bound = bounds(idx_bound);
        LL(idx_drift_scale, idx_bound) = log_joint_likelihood_3D([drift_scale, true_lambda, bound], drifts_cell, true_rt, fixed);
    end
end
figure
imagesc(LL)
xticks(1:res)
yticks(1:res)

xticklabels(bounds)
yticklabels(drifts)
%%

true_bound = 2.1;
res = 40;
drifts = linspace(1.0,1.6, res);
leaks = linspace(-0.1, 0.5, res);
LL = nan(res,res);
for idx_drift_scale = 1:res
    for idx_lambda = 1:res
        fprintf('done 1 \n')
        drift_scale = drifts(idx_drift_scale);
        lambda = leaks(idx_lambda);
        LL(idx_drift_scale, idx_lambda) = log_joint_likelihood_3D([drift_scale,lambda,true_bound], drifts_cell, true_rt, fixed);
    end
end
figure
imagesc(LL)
xticks(1:res)
yticks(1:res)

xticklabels(drifts)
yticklabels(leaks)



%% --- HELPER: Your Drift Generation Code ---
function drifts = generate_complex_drift(time_vec, dt, mean_drift)
time_steps = length(time_vec);

% Parameters (Hardcoded or passed as needed)
std_drift = 0.15; %0.2
correlation_time_ms = 100;
peak_amplitude = 6.5;
peak_probability = 0.005; %.001 % Slight adjust for density
peak_decay_ms = 20;
sampling_interval_ms = 20;
mean_std = .1;
peak_std = .1;

% 1. Correlated Noise
filter_length = round(correlation_time_ms / 1000 / dt);
if filter_length < 1, filter_length = 1; end
filter_coeffs = ones(1, filter_length) / filter_length;

raw = randn(1, time_steps);
correlated_noise = filter(filter_coeffs, 1, raw);
% Normalize std
correlated_noise = (correlated_noise - mean(correlated_noise)) / std(correlated_noise);

base_process = mean_drift + randn(1)*mean_std + std_drift*correlated_noise;

% 2. Peaks
decay_constant = round(peak_decay_ms / 1000 / dt);
if decay_constant < 1, decay_constant = 1; end

peak_triggers = find(rand(1, time_steps) < peak_probability);
peak_train = zeros(1, time_steps);

for trigger_idx = peak_triggers
    decay_length = min(5*decay_constant, time_steps - trigger_idx + 1);
    decay_shape = (1+randn(1)*peak_std)*exp(-(0:decay_length-1) / decay_constant);

    idx_end = trigger_idx + decay_length - 1;
    peak_train(trigger_idx : idx_end) = peak_train(trigger_idx : idx_end) + decay_shape;
end

if ~isempty(peak_triggers) && max(peak_train) > 0
    peak_train = peak_train / max(peak_train) * peak_amplitude;
end

drifts_continuous = base_process + peak_train;

% 3. Sample and Hold (Simulates discrete sampling)
sampling_interval_steps = round(sampling_interval_ms / 1000 / dt);
if sampling_interval_steps < 1, sampling_interval_steps = 1; end

drifts = zeros(1, time_steps);
for s_idx = 1 : sampling_interval_steps : time_steps
    measured = drifts_continuous(s_idx);
    e_idx = min(s_idx + sampling_interval_steps - 1, time_steps);
    drifts(s_idx : e_idx) = measured;
end

% Return as column vector for consistency
drifts = drifts(:);
end

%% --- OBJECTIVE FUNCTION ---
function ll = log_joint_likelihood_3D(params, drifts_cell, rts, fixed)
prop_mean_drift = params(1);
prop_lambda     = params(2);
prop_bound      = params(3); % NEW

% --- Dynamic Grid Setup ---
% We must calculate the grid parameters for THIS specific bound.
% x_min is fixed, dx is fixed.

% Calculate number of spatial nodes needed to reach this bound
prop_grid_size = round((prop_bound - fixed.x_min) / fixed.dx) + 1;

% Calculate where x0 (0.0) falls on this grid
prop_start_idx = round((fixed.x0 - fixed.x_min) / fixed.dx) + 1;

% Pack MEX params
mex_p = [fixed.dt, fixed.dx, fixed.sigma_sq, fixed.x0, ...
    fixed.x_min, prop_grid_size, prop_start_idx, prop_lambda];

N = length(rts);
log_liks = zeros(N, 1);

t_vec = (0:fixed.dt:fixed.T_max)';

idx_trunc = round(fixed.T_trunc / fixed.dt);

parfor i = 1:N
    current_drift = drifts_cell{i} * prop_mean_drift;

    % 1. Solve PDE (Get PDF and Final Survival Prob)
    [p_dist, survival_at_Tmax] = leaky_pde_robust(current_drift, mex_p);

    % 2. Calculate Normalization Factor (The Truncation Logic)
    % We need the total probability mass that crossed BEFORE T_trunc.
    % This is the area under the curve from 0 to T_trunc.
    if idx_trunc > 1
        % Sum flux up to truncation point
        prob_missing = sum(p_dist(1:idx_trunc)) * fixed.dt;
    else
        prob_missing = 0;
    end

    % The total probability of a valid observation in the truncated window
    normalization_factor = 1.0 - prob_missing;

    % Safety: If the model predicts EVERYTHING happens before T_trunc,
    % the normalization approaches 0.
    if normalization_factor < 1e-12
        normalization_factor = 1e-12;
    end

    % 3. Handle Observed Data

    % A. Censored Data (Subject didn't respond by Tmax)
    if rts(i) >= (fixed.T_max - fixed.dt) || isnan(rts(i))
        % Likelihood is Survival Probability renormalized
        lik = survival_at_Tmax / normalization_factor;

        % B. Valid Response Data
    else
        idx = round(rts(i) / fixed.dt);

        % Sanity check: If data violates truncation (impossible data)
        if idx < idx_trunc
            % This data point exists but is below the truncation threshold
            % defined in the likelihood. This implies a mismatch definition.
            lik = 1e-12;
        else
            % Standard lookup
            if idx > length(p_dist), idx = length(p_dist); end

            %                 raw_pdf = p_dist(idx);
            raw_pdf = interp1(t_vec,p_dist,rts(i));

            % Renormalize PDF
            lik = raw_pdf / normalization_factor;
        end
    end

    % 4. Numerical Stability
    if lik < 1e-12 || isnan(lik)
        lik = 1e-12;
    end

    log_liks(i) = log(lik);
end

ll = sum(log_liks);
end

%% --- PRIORS ---
function lprior = log_trapezoidal_prior(params, LB, UB, PLB, PUB)
% params: [1 x D] vector of parameters
% LB, UB: Hard bounds
% PLB, PUB: Plausible (Plateau) bounds

D = length(params);
lprior = 0;

penalty = 20.0;

for i = 1:D
    x = params(i);

    % Check Hard Bounds (Immediate rejection)
    if x < LB(i) || x > UB(i)
        lprior = -inf;
        return;
    end

    % Calculate Trapezoidal shape
    if x >= PLB(i) && x <= PUB(i)
        % On the plateau
        val = 0;
    elseif x < PLB(i)
        % Left Slope (Linear decay from PLB to LB)
        % We map [LB, PLB] -> [-Inf, 0] ?
        % Actually, standard trapezoidal usually implies linear probability density.
        % Linear PDF = Linear Slope. Log PDF = Log(Linear).
        % However, usually in optimization, we prefer "Smoothed Box" priors.

        % Let's use a Linear Log-Prior (Exponential decay in probability)
        % This acts like a soft constraint.
        % slope * (distance from plateau)

        % Normalized distance from plateau (0 at PLB, 1 at LB)
        dist = (PLB(i) - x) / (PLB(i) - LB(i));

        % Penalty factor (arbitrary, but ensures -10 or -20 at the edge)
        val = -penalty * dist;

    else % x > PUB(i)
        % Right Slope
        dist = (x - PUB(i)) / (UB(i) - PUB(i));
        val = -penalty * dist;
    end

    lprior = lprior + val;
end
end

function lprior = log_correlated_prior_3D(params, LB, UB, PLB, PUB)
% Params: [Mean_Drift, Lambda, Bound]
v      = params(1);
lambda = params(2);
bound  = params(3);

% --- 1. Hard Bounds Check ---
if any(params < LB) || any(params > UB)
    lprior = -inf;
    return;
end

% --- 2. Independent Priors for Drift and Bound ---
% We use the Trapezoidal logic for these two, as they are the "independent" drivers.

% Prior on Drift (Index 1)
lp_v = log_trapezoidal_1D(v, LB(1), UB(1), PLB(1), PUB(1));

% Prior on Bound (Index 3)
lp_b = log_trapezoidal_1D(bound, LB(3), UB(3), PLB(3), PUB(3));

% --- 3. Conditional Prior on Lambda ---
% Physics Constraint: Asymptote (v/lambda) should be > Bound.
% Let's target an asymptote of roughly 1.3 * Bound.
% Therefore: Target_Lambda = v / (1.3 * Bound).

target_ratio = 1.3;

% Avoid division by zero if bound is crazy small (unlikely due to LB)
mu_lambda = v / (target_ratio * bound);

% How strict?
% Sigma = 0.5 allows the asymptote to vary roughly between 1.0*B and 2.0*B
sigma_lambda = 0.75;

% Log-Normal-ish penalty (Gaussian in linear space centered on mu)
lp_lambda = -0.5 * ((lambda - mu_lambda) / sigma_lambda)^2;

% --- 4. Total ---
lprior = lp_v + lp_b + lp_lambda;
end

function lprior = log_correlated_prior_3D_robust(params, LB, UB, PLB, PUB)
% Params: [Drift, Lambda, Bound]
v      = params(1);
lambda = params(2);
bound  = params(3);

% 1. Hard Bounds Check (Immediate exit)
if any(params < LB) || any(params > UB)
    lprior = -inf; return;
end

% 2. Independent Priors for Drift and Bound (Trapezoidal)
lp_v = log_trapezoidal_1D(v, LB(1), UB(1), PLB(1), PUB(1));
lp_b = log_trapezoidal_1D(bound, LB(3), UB(3), PLB(3), PUB(3));

% 3. The New Conditional Prior on Lambda
% We define a "Safe Limit" for Lambda.
% If lambda is below this limit, the prior is flat (0 penalty).
% If lambda exceeds this limit, we penalize.

% Constraint: Asymptote should be at least 1.1x the Bound.
% v / lambda > 1.1 * bound
% lambda < v / (1.1 * bound)

safety_margin = 1.1;

% Avoid div by zero if v is tiny (though usually v > 0)
% If v=0, lambda must be 0.
if v < 1e-3
    lambda_limit = 0;
else
    lambda_limit = v / (safety_margin * bound);
end

if lambda <= lambda_limit
    % SAFE ZONE: Lambda is small enough that the particle
    % easily reaches the bound.
    % We add a weak preference for smaller lambdas to break degeneracy at 0?
    % No, let's keep it flat to allow perfect integration (lambda=0) unbiased.
    lp_lambda = 0;

else
    % DANGER ZONE: Lambda is too high. Asymptote is below the bound.
    % Soft penalty (half-Gaussian)
    sigma_penalty = 0.75;
    lp_lambda = -0.5 * ((lambda - lambda_limit) / sigma_penalty)^2;
end

% 4. Total
lprior = lp_v + lp_b + lp_lambda;
end

% Helper for single trapezoid
function val = log_trapezoidal_1D(x, lb, ub, plb, pub)
if x >= plb && x <= pub
    val = 0;
elseif x < plb
    % Normalized distance to left plateau edge
    dist = (plb - x) / (plb - lb);
    val = -10.0 * dist; % Penalty of -10 at the hard edge
else % x > pub
    dist = (x - pub) / (ub - pub);
    val = -10.0 * dist;
end
end