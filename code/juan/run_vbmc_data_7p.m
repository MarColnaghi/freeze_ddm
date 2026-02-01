function run_vbmc_data_7p()
% this version includes two separate betas for the drift, one inside the
% loom evoked window and one outside

    % Check for VBMC
    if isempty(which('vbmc'))
        error('VBMC is not in the path. Please install/add it.');
    end

    %% 1. Setup Parameters
    rng(42); 
    
    % increased VBMC accuracy at the expense of time?
    inc_vbmc_acc = 0;
   
    
    % -- Likelihood / Signal Parameters (Coarse) --
    fixed.dt = 1/300;       % 4ms for signal and PDE
    fixed.dx = 0.025;
    fixed.sigma_sq = 1.0;
%     fixed.bound = 1.0;
    fixed.x0 = 0.0;
    
    % -- Simulation Parameters (Fine) --
    sim.dt = 0.00005;       % 0.05ms for ground truth simulation
    
    % PDE Grid Setup
    x_min = -7; 
    
    % Save for objective function
    fixed.x_min = x_min;
    
    % Truncation and T_max
    fixed.T_trunc = .3;
    fixed.T_trunc = .3;
    fixed.T_max = 10.5; 
    
    n_moving = 0:4;
    sloom = 25;
    nloom_min = 8;
    nloom_max = 20;

    %% 2. Data processing
    
    tic;
    
    % load data
    load("G:\My Drive\MATLAB CODE\FlyFreeze\Code_from_Marco\soc_mot_array.mat");
    load("G:\My Drive\MATLAB CODE\FlyFreeze\Code_from_Marco\bouts_proc.mat");
    
    % select data
    idx = ismember(bouts_proc.moving_flies,n_moving) & bouts_proc.sloom == sloom ...
        & bouts_proc.durations_s >= fixed.T_trunc ...
        & bouts_proc.nloom>=nloom_min ...
        & bouts_proc.nloom<=nloom_max;
    median_fs = median(bouts_proc.avg_fs_1s(idx));
    % filter by focal speed
    idx = idx & (bouts_proc.avg_fs_1s < median_fs);
    soc_mot_array = soc_mot_array(idx,:);
    bouts_proc = bouts_proc(idx,:);
    
    N_trials = size(soc_mot_array,1);
    
    sim.dt = 1/60;        
        
    % Time vectors
    t_coarse = 0:sim.dt:fixed.T_max;
    t_fine   = 0:fixed.dt:fixed.T_max;
    
    fprintf('\n %d trials loaded \n', N_trials);
    
    
%     fine_signal; %initialize
    
    drifts_cell = cell(N_trials, 1); % Stores the COARSE signal (what the model sees)
    true_rt = bouts_proc.durations_s(1:N_trials);
        
    for i = 1:N_trials
        
        sim_signal_coarse = [soc_mot_array(i,:),soc_mot_array(i,end)];
        drifts_cell{i} = interp1(t_coarse(:), sim_signal_coarse(:), t_fine(:), 'nearest');
        
    end

    fprintf('Data processing complete in %.2f seconds.\n', toc);

    %% 3. Setup VBMC
    
    % Order: [Drift, Lambda, Bound]
    % Note: Bound should be positive.
    
    %       Drift1  Drift2   Lambda   Bound  ndt  lambda  pcont
    LB  = [ 0.01,   0.01,    -.5,     0.3    0      .1    0];
    UB  = [ 6.0,    6.0,     3.0,     4.0    .5      50   .1];
    
    % Plausible
    PLB = [ 0.25,   0.25,    -.05,    0.5    .005    .2    .01];
    PUB = [ 5.0,    5.0,     2.0,     2.5    .3       20    .09];
    
    x0_start = [2, 1, .1, 1, .05, 4, .05]; 

    % Update Handles
%     log_prior_fun = @(p) log_correlated_prior_3D_robust(p, LB, UB, PLB, PUB);
    log_prior_fun = @(p) log_trapezoidal_prior(p, LB, UB, PLB, PUB);
    
    % Note: Likelihood function now needs to handle the dynamic bound internally
    log_lik_fun = @(p) log_joint_likelihood_7D(p, drifts_cell, true_rt, fixed);

    % C. Define The Posterior (Target for VBMC)
    % Posterior = Likelihood + Prior
    log_posterior_fun = @(p) log_lik_fun(p) + log_prior_fun(p);

    options = vbmc('defaults');
    options.Display = 'iter';
    options.Plot = false; 
    options.SpecifyTargetNoise = false;
    options.TolStableCount = 80;
    options.MinFinalComponents = 50;
    options.MaxFunEvals = 500;
    
    % EXTRA OPTIONS FOR INCREASED ACCURACY
    if inc_vbmc_acc
        % 1. Force a high-resolution initial map
        options.FunEvalStart = 50;
        % 2. Force it to run longer (don't stop on stability early)
        options.MinFunEvals = 400;
        options.MaxFunEvals = 800;
        % 3. Force higher precision before declaring convergence
        options.TolImprovement = 0.002;
        % 4. Start with a more complex mixture model (smoother tails)
        options.Kwarmup = 5;
    end
    

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
    fprintf('Drift 1:     %.2f (+/- %.2f)\n', post_mean(1), post_std(1));
    fprintf('Drift 2:     %.2f (+/- %.2f)\n', post_mean(2), post_std(2));
    fprintf('Lambda:      %.2f (+/- %.2f)\n', post_mean(3), post_std(3));
    fprintf('Bound:       %.2f (+/- %.2f)\n', post_mean(4), post_std(4));
    fprintf('Ndt:         %.2f (+/- %.2f)\n', post_mean(5), post_std(5));
    fprintf('cont. lam.:  %.2f (+/- %.2f)\n', post_mean(6), post_std(6));
    fprintf('cont. prob.: %.2f (+/- %.2f)\n', post_mean(7), post_std(7));
    
    vbmc_plot(vp); sgtitle('Parameter Recovery');
    
    %obtained log likelihood
    obtained_LL = log_joint_likelihood_7D(post_mean, drifts_cell, true_rt, fixed);
    %true log likelihood
%     true_LL = log_joint_likelihood_3D([true_drift_scale,true_lambda,true_bound], drifts_cell, true_rt, fixed);
    fprintf('Obtained Log Likelihood: %.2f\n', obtained_LL);
%     fprintf('True Log Likelihood:     %.2f\n', true_LL);
   
    
    
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

%%%%%%%%%%%%%%%%
% with ndt
function ll = log_joint_likelihood_7D(params, drifts_cell, rts, fixed)
    % Unpack
    drift1          = params(1);
    drift2          = params(2);
    prop_lambda     = params(3);
    prop_bound      = params(4);
    prop_tnd        = params(5);
    lambda_cont     = params(6); % Rate of exponential (1/mean)
    p_cont          = params(7); % Mixture probability
    
    % --- Dynamic Grid Setup ---
    prop_grid_size = round((prop_bound - fixed.x_min) / fixed.dx) + 1;
    prop_start_idx = round((fixed.x0 - fixed.x_min) / fixed.dx) + 1;
    
    % Pack MEX params
    mex_p = [fixed.dt, fixed.dx, fixed.sigma_sq, fixed.x0, ...
             fixed.x_min, prop_grid_size, prop_start_idx, prop_lambda];
         
    N = length(rts);
    log_liks = zeros(N, 1);
    
    % Time vector for DDM decision times (0 to T_max)
    t_pde = (0 : length(drifts_cell{1})-1)' * fixed.dt;
    
    % Pre-calculate Contaminant Normalization
    % Contaminants are defined on t > 0, but we truncate observations < T_trunc.
    % Mass of Exponential CDF at T_trunc:
    if fixed.T_trunc > 0
        missing_mass_cont = 1 - exp(-lambda_cont * fixed.T_trunc);
    else
        missing_mass_cont = 0;
    end
    
    idx_05s = floor(0.5 / fixed.dt);
    
    parfor i = 1:N
        current_drift = drifts_cell{i};
        current_drift(1:idx_05s) = current_drift(1:idx_05s) * drift1;
        current_drift((idx_05s+1):end) = current_drift((idx_05s+1):end) * drift2;
        
        % 1. Solve PDE (DDM Component Only)
        [p_dist, ~] = leaky_pde_robust(current_drift, mex_p);
        
        % --- 2. Calculate DDM Normalization ---
        % DDM "starts" at T_nd.
        % We miss any DDM mass that crosses before (T_trunc - T_nd).
        
        trunc_decision_time = fixed.T_trunc - prop_tnd;
        
        if trunc_decision_time <= 0
            % T_nd is large; DDM starts after the truncation window.
            missing_mass_ddm = 0;
        else
            idx_trunc = floor(trunc_decision_time / fixed.dt);
            if idx_trunc < 1
                missing_mass_ddm = 0;
            elseif idx_trunc >= length(p_dist)
                missing_mass_ddm = 1.0;
            else
                missing_mass_ddm = sum(p_dist(1:idx_trunc)) * fixed.dt;
            end
        end
        
        % COMBINED Normalization Factor
        % Z = (1-pi)*Z_ddm + pi*Z_cont
        term_ddm  = (1 - p_cont) * (1 - missing_mass_ddm);
        term_cont = p_cont       * (1 - missing_mass_cont);
        
        normalization_factor = term_ddm + term_cont;
        
        if normalization_factor < 1e-12, normalization_factor = 1e-12; end

        % --- 3. Handle Observed Data ---
        
        % A. Censored Data (RT >= T_max or NaN)
        if rts(i) >= (fixed.T_max - fixed.dt) || isnan(rts(i))
            
            % 1. DDM Survival: Prob(DecisionTime > T_max - T_nd)
            cutoff_decision = fixed.T_max - prop_tnd;
            
            if cutoff_decision <= 0
                surv_ddm = 1.0; % DDM hasn't started yet
            else
                idx_cut = floor(cutoff_decision / fixed.dt);
                if idx_cut >= length(p_dist)
                    absorbed = 1.0;
                elseif idx_cut < 1
                    absorbed = 0.0;
                else
                    absorbed = sum(p_dist(1:idx_cut)) * fixed.dt;
                end
                surv_ddm = 1.0 - absorbed;
            end
            
            % 2. Contaminant Survival: Prob(T > T_max) = exp(-lambda * T_max)
            surv_cont = exp(-lambda_cont * fixed.T_max);
            
            % Mixture
            lik = ((1 - p_cont) * surv_ddm + p_cont * surv_cont) / normalization_factor;
            
        % B. Valid Response Data
        else
            % 1. DDM Likelihood
            % Evaluated at (RT - T_nd)
            decision_time = rts(i) - prop_tnd;
            
            if decision_time <= 0
                % RT is faster than T_nd. DDM prob is strictly 0.
                lik_ddm = 0;
            else
                % Lookup
                raw_ddm = interp1(t_pde, p_dist, decision_time, 'linear', 0);
                
                % Check if valid wrt truncation
                if decision_time < trunc_decision_time
                    lik_ddm = 0;
                else
                    lik_ddm = raw_ddm;
                end
            end
            
            % 2. Contaminant Likelihood
            % Evaluated at RT (Not shifted!)
            % PDF = lambda * exp(-lambda * t)
            lik_cont = lambda_cont * exp(-lambda_cont * rts(i));
            
            % Mixture
            raw_mixture = (1 - p_cont) * lik_ddm + p_cont * lik_cont;
            lik = raw_mixture / normalization_factor;
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