function plot_fit_10p()


    %% 1. Setup Parameters
    rng(42); 
   
    
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
    fixed.T_trunc = .5;
    fixed.T_max = 10.5; 
    
    % fit params
    % fit_params = [3.40, 2.01, 1.63, .59, .11, .60, .01, .44, 8.71, .06];
    % fit_params = [1.2, -.19, .51, 3.56, 2.42, 1.57, .00, .59, .28, .10];
    % fit_params = [1.17, -.18, .36, 3.47, 3.00, 1.40, .01, .64, .32, .09];
    % fit_params = [3.04, 3.97, 1.09, 0.87, 0.06, 0.24, .00, .17, 4.48, .03];
    % fit_params = [3.46, 4.95, 1.02, 0.97, 0.05, 0.26, .00, .11, 4.52, .03];
    % fit_params = [3.30, 5.03, 1.07, 0.69, 0.06, 0.82, .00, .56, 4.54, .03];
    fit_params = [3.88, 0, 1.14, 0.74, 0.08, 0.58, .00, .35, 4.51, .03]; %4.81


    % fit_params = [2.10, 1.02, 0.07, 1.12, .01, 6.77, 0.08];
    n_moving = 0:4;
    sloom = 50;
    nloom_min = 2;
    nloom_max = 20;
    

    %% 2. Data processing
    
    tic;
    
    % load data
    load("G:\My Drive\MATLAB CODE\FlyFreeze\Code_from_Marco\soc_mot_array.mat");
    load("G:\My Drive\MATLAB CODE\FlyFreeze\Code_from_Marco\bouts_proc.mat");
    
    % select data
    idx = ismember(bouts_proc.moving_flies,n_moving) & bouts_proc.sloom == sloom ...
         & bouts_proc.durations_s >= fixed.T_trunc...
         & bouts_proc.nloom>=nloom_min ...
         & bouts_proc.nloom<=nloom_max;
    median_fs = median(bouts_proc.avg_fs_1s(idx));
    % filter by focal speed
    idx = idx & (bouts_proc.avg_fs_1s >= median_fs);
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
        
        sim_signal_coarse = [0,soc_mot_array(i,:)];
        drifts_cell{i} = interp1(t_coarse(:), sim_signal_coarse(:), t_fine(:), 'nearest');
        
    end

    fprintf('Data processing complete in %.2f seconds.\n', toc);

    %% 3. Plot histogram
    
    figure('Position', [100, 100, 900, 600], 'color', 'w');
    hold on;
    t_vec = fixed.T_trunc:(5/60):20;
    % Plot the histogram of simulated FPTs
    % 'Normalization', 'pdf' ensures the total area of the histogram is 1.
    
    histogram(true_rt,t_vec,'Normalization', 'pdf', ...
        'DisplayName', 'Data',...
        'LineWidth',1.5,'DisplayStyle', 'stairs');
    
    xlim([0,fixed.T_max]);
   
      
    % combined likelihood
    
    combined_lik = combined_likelihood_10D(fit_params, drifts_cell, fixed);
%     t_vec = (0:fixed.dt:fixed.T_max)';
    t_vec = (0:fixed.dt:fixed.T_max)';
    
    plot(t_vec,combined_lik, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Fit');
    xlabel('Freeze duration (s)'); ylabel('Probability Density');
    legend('show', 'Location', 'northeast'); grid on; box on; axis tight;
    
    obtained_LL = log_joint_likelihood_10D(fit_params, drifts_cell, true_rt, fixed);
    fprintf('Obtained Log Likelihood: %.2f\n', obtained_LL);
    
    % fit_params = [1.62, 1.29, 0.04, .62, .01, 1, 0.05];
    fit_params = [3.39, 1.97, 2.04, 1.30, .01, 2.17, 0.09];
    combined_lik3 = combined_likelihood_7D(fit_params, drifts_cell, fixed);   
    % plot(t_vec,combined_lik3, 'k-', 'LineWidth', 1, 'DisplayName', 'Fit');
    
    xlim([0,fixed.T_max]);
    
    obtained_LL = log_joint_likelihood_7D(fit_params, drifts_cell, true_rt, fixed);
    fprintf('Obtained Log Likelihood: %.2f\n', obtained_LL);
    
end





function lik = combined_likelihood_7D(params, drifts_cell, fixed)
    % Unpack
    drift1          = params(1);
    drift2          = params(2);
    prop_lambda     = params(3);
    prop_bound      = params(4);
    prop_tnd        = params(5); % NEW: Non-Decision Time
    lambda_cont     = params(6); % Rate of exponential (1/mean)
    p_cont          = params(7); % Mixture probability
    

    % --- Dynamic Grid Setup ---
    prop_grid_size = round((prop_bound - fixed.x_min) / fixed.dx) + 1;
    prop_start_idx = round((fixed.x0 - fixed.x_min) / fixed.dx) + 1;
    
    % Pack MEX params
    mex_p = [fixed.dt, fixed.dx, fixed.sigma_sq, fixed.x0, ...
             fixed.x_min, prop_grid_size, prop_start_idx, prop_lambda];
         
    N = length(drifts_cell);    
    t_vec = (0:fixed.dt:fixed.T_max)';
    
    idx_trunc = round(fixed.T_trunc / fixed.dt);
    
    lik = zeros(size(t_vec));
    
    % Pre-calculate Contaminant Normalization
    % Contaminants are defined on t > 0, but we truncate observations < T_trunc.
    % Mass of Exponential CDF at T_trunc:
    if fixed.T_trunc > 0
        missing_mass_cont = 1 - exp(-lambda_cont * fixed.T_trunc);
    else
        missing_mass_cont = 0;
    end
    
    % contaminant
    lik_cont = lambda_cont * exp(-lambda_cont * t_vec);
    lik_cont(1:idx_trunc) = 0;
    
    idx_05s = floor(0.6 / fixed.dt);
    
    for i = 1:N
        current_drift = drifts_cell{i};
        current_drift(1:idx_05s) = current_drift(1:idx_05s) * drift1;
        current_drift((idx_05s+1):end) = current_drift((idx_05s+1):end) * drift2;
        
        % 1. Solve PDE
        % p_dist represents the distribution of DECISION times
        % survival_at_end is mass remaining after the full drift vector duration
        [p_ddm, ~] = leaky_pde_robust(current_drift, mex_p);
        
        
        % --- 2. Calculate Normalization Factor (Truncation) ---
        % We only observe RTs > T_trunc.
        % This implies we only observe Decision Times > (T_trunc - T_nd).
        
        trunc_decision_time = fixed.T_trunc - prop_tnd;
        
        idx_trunc = 1;
        if trunc_decision_time <= 0
            % Case: T_nd is larger than T_trunc. 
            % The diffusion process effectively starts "after" the truncation window.
            % No part of the PDE distribution is truncated.
            missing_mass_ddm = 0;
        else
            % Case: We need to remove the probability mass that occurred 
            % before the effective truncation time.
            idx_trunc = floor(trunc_decision_time / fixed.dt);
            
            if idx_trunc < 1
                missing_mass_ddm = 0;
            elseif idx_trunc >= length(p_ddm)
                missing_mass_ddm = 1.0;
            else
                missing_mass_ddm = sum(p_ddm(1:idx_trunc)) * fixed.dt;
            end
                 
        end
        
        % COMBINED Normalization Factor
        % Z = (1-pi)*Z_ddm + pi*Z_cont
        term_ddm  = (1 - p_cont) * (1 - missing_mass_ddm);
        term_cont = p_cont       * (1 - missing_mass_cont);
        
        normalization_factor = term_ddm + term_cont;
        
        if normalization_factor < 1e-12, normalization_factor = 1e-12; end
        
        idx_tnd = floor(prop_tnd/fixed.dt);
        
        % truncate
        p_ddm(idx_tnd+1:end) = p_ddm(1:length(t_vec)-idx_tnd);
        p_ddm(1:idx_tnd) = 0;        
        p_ddm(1:idx_trunc) = 0;
        
        
            
        % Mixture
        raw_mixture = (1 - p_cont) * p_ddm + p_cont * lik_cont;
        
        % normalize
        p_dist = raw_mixture / normalization_factor;
        
        lik = lik + p_dist/N;
        
        lik(lik < 1e-12 | isnan(lik)) = 1e-12;

    end
end


function lik = combined_likelihood_10D(params, drifts_cell, fixed)
    % Unpack
    drift1          = params(1);
    prop_lambda1    = params(2);
    prop_bound1     = params(3);
    drift2          = params(4);
    prop_lambda2    = params(5);
    prop_bound2     = params(6);
    prop_tnd        = params(7);
    pmix            = params(8);
    lambda_cont     = params(9); % Rate of exponential (1/mean)
    p_cont          = params(10); % Mixture probability
    
    % --- Dynamic Grid Setup ---
    prop_grid_size1 = round((prop_bound1 - fixed.x_min) / fixed.dx) + 1;
    prop_grid_size2 = round((prop_bound2 - fixed.x_min) / fixed.dx) + 1;
    prop_start_idx = round((fixed.x0 - fixed.x_min) / fixed.dx) + 1;
    
    % Pack MEX params
    mex_p1 = [fixed.dt, fixed.dx, fixed.sigma_sq, fixed.x0, ...
             fixed.x_min, prop_grid_size1, prop_start_idx, prop_lambda1];
    mex_p2 = [fixed.dt, fixed.dx, fixed.sigma_sq, fixed.x0, ...
             fixed.x_min, prop_grid_size2, prop_start_idx, prop_lambda2];
         
    N = length(drifts_cell);   
    
    idx_trunc = round(fixed.T_trunc / fixed.dt);
    
    % Time vector for DDM decision times (0 to T_max)
    t_vec = (0:fixed.dt:fixed.T_max)';

    lik = zeros(size(t_vec));
    
    % Pre-calculate Contaminant Normalization
    % Contaminants are defined on t > 0, but we truncate observations < T_trunc.
    % Mass of Exponential CDF at T_trunc:
    if fixed.T_trunc > 0
        missing_mass_cont = 1 - exp(-lambda_cont * fixed.T_trunc);
    else
        missing_mass_cont = 0;
    end
    
    % contaminant
    lik_cont = lambda_cont * exp(-lambda_cont * t_vec);
    lik_cont(1:idx_trunc) = 0;
    
    % idx_05s = floor(0.5 / fixed.dt);
    
    for i = 1:N
        current_drift1 = drifts_cell{i} * drift1;
        current_drift2 = drifts_cell{i} * drift2;
        
        
        % 1. Solve PDE (DDM Component Only)
        [p_dist1, ~] = leaky_pde_robust(current_drift1, mex_p1);
        [p_dist2, ~] = leaky_pde_robust(current_drift2, mex_p2);
        
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
            elseif idx_trunc >= length(p_dist1)
                missing_mass_ddm = 1.0;
            else
                missing_mass_ddm = pmix * sum(p_dist1(1:idx_trunc)) * fixed.dt +...
                                   (1-pmix) * sum(p_dist2(1:idx_trunc)) * fixed.dt;
            end
        end
        
        % COMBINED Normalization Factor
        % Z = (1-pi)*Z_ddm + pi*Z_cont
        term_ddm  = (1 - p_cont) * (1 - missing_mass_ddm);
        term_cont = p_cont       * (1 - missing_mass_cont);
        
        normalization_factor = term_ddm + term_cont;
        
        if normalization_factor < 1e-12, normalization_factor = 1e-12; end

        
        idx_tnd = floor(prop_tnd/fixed.dt);

        p_ddm = pmix*p_dist1 + (1-pmix)*p_dist2;
        
        % truncate
        p_ddm(idx_tnd+1:end) = p_ddm(1:length(t_vec)-idx_tnd);
        p_ddm(1:idx_tnd) = 0;        
        p_ddm(1:idx_trunc) = 0;
                          
        % Mixture
        raw_mixture = (1 - p_cont) * p_ddm + p_cont * lik_cont;
        
        % normalize
        p_dist = raw_mixture / normalization_factor;
        
        lik = lik + p_dist/N;
        
        lik(lik < 1e-12 | isnan(lik)) = 1e-12;

    end

end


function ll = log_joint_likelihood_10D(params, drifts_cell, rts, fixed)
    % Unpack
    drift1          = params(1);
    prop_lambda1    = params(2);
    prop_bound1     = params(3);
    drift2          = params(4);
    prop_lambda2    = params(5);
    prop_bound2     = params(6);
    prop_tnd        = params(7);
    pmix            = params(8);
    lambda_cont     = params(9); % Rate of exponential (1/mean)
    p_cont          = params(10); % Mixture probability
    
    % --- Dynamic Grid Setup ---
    prop_grid_size1 = round((prop_bound1 - fixed.x_min) / fixed.dx) + 1;
    prop_grid_size2 = round((prop_bound2 - fixed.x_min) / fixed.dx) + 1;
    prop_start_idx = round((fixed.x0 - fixed.x_min) / fixed.dx) + 1;
    
    % Pack MEX params
    mex_p1 = [fixed.dt, fixed.dx, fixed.sigma_sq, fixed.x0, ...
             fixed.x_min, prop_grid_size1, prop_start_idx, prop_lambda1];
    mex_p2 = [fixed.dt, fixed.dx, fixed.sigma_sq, fixed.x0, ...
             fixed.x_min, prop_grid_size2, prop_start_idx, prop_lambda2];
         
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
    
    % idx_05s = floor(0.5 / fixed.dt);
    
    for i = 1:N
        current_drift1 = drifts_cell{i} * drift1;
        current_drift2 = drifts_cell{i} * drift2;
        
        
        % 1. Solve PDE (DDM Component Only)
        [p_dist1, ~] = leaky_pde_robust(current_drift1, mex_p1);
        [p_dist2, ~] = leaky_pde_robust(current_drift2, mex_p2);
        
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
            elseif idx_trunc >= length(p_dist1)
                missing_mass_ddm = 1.0;
            else
                missing_mass_ddm = pmix * sum(p_dist1(1:idx_trunc)) * fixed.dt +...
                                   (1-pmix) * sum(p_dist2(1:idx_trunc)) * fixed.dt;
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
                if idx_cut >= length(p_dist1)
                    absorbed = 1.0;
                elseif idx_cut < 1
                    absorbed = 0.0;
                else
                    absorbed = pmix * sum(p_dist1(1:idx_cut)) * fixed.dt +...
                               (1-pmix) * sum(p_dist2(1:idx_cut)) * fixed.dt;
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
                raw_ddm1 = interp1(t_pde, p_dist1, decision_time, 'linear', 0);
                raw_ddm2 = interp1(t_pde, p_dist2, decision_time, 'linear', 0);
                
                % Check if valid wrt truncation
                if decision_time < trunc_decision_time
                    lik_ddm = 0;
                else
                    lik_ddm = pmix*raw_ddm1 + (1-pmix)*raw_ddm2;
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
    
    idx_05s = floor(0.6 / fixed.dt);
    
    for i = 1:N
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