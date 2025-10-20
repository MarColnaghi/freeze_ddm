function lik = extrema_detection_likelihood(signal, rt, bound, dt, t_nd_mean, t_nd_std)
    % Calculate the likelihood of observing a reaction time under the extrema detection model.
    %
    % Parameters:
    % -----------
    % signal : vector (n_samples x 1)
    %     Time-varying signal for a single trial. Should be discretely sampled.
    % rt : scalar
    %     Observed reaction time (in seconds)
    % bound : scalar
    %     Detection threshold (positive value)
    % dt : scalar (optional)
    %     Time step between samples (in seconds). Default 0.001 (1ms)
    % t_nd_mean : scalar (optional)
    %     Mean of non-decision time (in seconds). Default 0.3
    % t_nd_std : scalar (optional)
    %     Standard deviation of non-decision time (in seconds). Default 0.05
    %
    % Returns:
    % --------
    % lik : scalar
    %     Probability density of observing this RT under the model
    %
    % Notes:
    % ------
    % The model assumes:
    % - An extremum is detected when signal(i) > bound
    % - Decision time is when first extremum detected
    % - RT = decision_time + non_decision_time
    % - Non-decision time ~ Normal(t_nd_mean, t_nd_std)
    
    % Set defaults
    if nargin < 4 || isempty(dt)
        dt = 0.001;
    end
    if nargin < 5 || isempty(t_nd_mean)
        t_nd_mean = 0.3;
    end
    if nargin < 6 || isempty(t_nd_std)
        t_nd_std = 0.05;
    end
    
    signal = signal(:);  % Ensure column vector
    n_samples = length(signal);
    
    % Calculate probability of exceeding bound at each time step
    % P(X > bound) where X ~ N(signal(i), sqrt(dt))
    p_exceed = normcdf(bound, signal, sqrt(dt), 'upper');
    
    % Probability of NOT exceeding bound up to each time point
    p_not_exceed = 1 - p_exceed;
    
    % Probability of first detection at each time step
    % P(detect at time t) = P(exceed at t) * P(didn't exceed before t)
    p_first_detection = zeros(n_samples, 1);
    cumulative_survival = 1.0;
    
    for i = 1:n_samples
        p_first_detection(i) = p_exceed(i) * cumulative_survival;
        cumulative_survival = cumulative_survival * p_not_exceed(i);
    end
    
    % Time grid (in seconds)
    time_grid = (0:n_samples-1)' * dt;
    
    % For each possible decision time, calculate likelihood of observed RT
    % RT = decision_time + non_decision_time
    % So: non_decision_time = RT - decision_time
    
    lik = 0.0;
    for i = 1:n_samples
        decision_time = time_grid(i);
        implied_nd_time = rt - decision_time;
        
        % Probability density of this non-decision time
        p_nd = normpdf(implied_nd_time, t_nd_mean, t_nd_std);
        
        % Total contribution to likelihood
        lik = lik + p_first_detection(i) * p_nd;
    end
end


function log_lik = log_likelihood_dataset(signals, rts, bound, dt, t_nd_mean, t_nd_std)
    % Calculate log-likelihood for entire dataset.
    %
    % Parameters:
    % -----------
    % signals : cell array of vectors
    %     Time-varying signal for each trial
    % rts : vector (n_trials x 1)
    %     Observed reaction times for each trial
    % bound : scalar
    %     Detection threshold
    % dt : scalar (optional)
    %     Time step between samples. Default 0.001
    % t_nd_mean : scalar (optional)
    %     Mean of non-decision time. Default 0.3
    % t_nd_std : scalar (optional)
    %     Standard deviation of non-decision time. Default 0.05
    %
    % Returns:
    % --------
    % log_lik : scalar
    %     Sum of log-likelihoods across all trials
    
    % Set defaults
    if nargin < 4 || isempty(dt)
        dt = 0.001;
    end
    if nargin < 5 || isempty(t_nd_mean)
        t_nd_mean = 0.3;
    end
    if nargin < 6 || isempty(t_nd_std)
        t_nd_std = 0.05;
    end
    
    n_trials = length(rts);
    log_lik = 0.0;
    
    for i = 1:n_trials
        lik = extrema_detection_likelihood(signals{i}, rts(i), bound, dt, t_nd_mean, t_nd_std);
        % Add small epsilon to avoid log(0)
        log_lik = log_lik + log(lik + 1e-300);
    end
end


function [params, neg_log_lik, exitflag] = fit_extrema_model(signals, rts, dt, varargin)
    % Fit extrema detection model parameters to data using maximum likelihood.
    %
    % Parameters:
    % -----------
    % signals : cell array of vectors
    %     Time-varying signal for each trial
    % rts : vector (n_trials x 1)
    %     Observed reaction times for each trial
    % dt : scalar
    %     Time step between samples (in seconds)
    % 
    % Optional name-value pairs:
    % ---------------------------
    % 'InitBound' : scalar (default: 0.1)
    %     Initial guess for detection bound
    % 'InitTndMean' : scalar (default: 0.3)
    %     Initial guess for mean non-decision time
    % 'InitTndStd' : scalar (default: 0.05)
    %     Initial guess for std of non-decision time
    % 'LowerBounds' : 3-element vector (default: [0.01, 0.1, 0.01])
    %     Lower bounds for [bound, t_nd_mean, t_nd_std]
    % 'UpperBounds' : 3-element vector (default: [1.0, 0.5, 0.2])
    %     Upper bounds for [bound, t_nd_mean, t_nd_std]
    %
    % Returns:
    % --------
    % params : struct with fields
    %     bound : fitted detection threshold
    %     t_nd_mean : fitted mean non-decision time
    %     t_nd_std : fitted std of non-decision time
    % neg_log_lik : scalar
    %     Negative log-likelihood at optimum
    % exitflag : integer
    %     Exit flag from fmincon
    
    % Parse optional inputs
    p = inputParser;
    addParameter(p, 'InitBound', 0.1);
    addParameter(p, 'InitTndMean', 0.3);
    addParameter(p, 'InitTndStd', 0.05);
    addParameter(p, 'LowerBounds', [0.01, 0.1, 0.01]);
    addParameter(p, 'UpperBounds', [1.0, 0.5, 0.2]);
    parse(p, varargin{:});
    
    % Initial parameter vector: [bound, t_nd_mean, t_nd_std]
    x0 = [p.Results.InitBound, p.Results.InitTndMean, p.Results.InitTndStd];
    lb = p.Results.LowerBounds;
    ub = p.Results.UpperBounds;
    
    % Objective function: negative log-likelihood
    obj_fun = @(x) -log_likelihood_dataset(signals, rts, x(1), dt, x(2), x(3));
    
    % Set up optimization options
    options = optimoptions('fmincon', ...
        'Display', 'iter', ...
        'Algorithm', 'interior-point', ...
        'MaxFunctionEvaluations', 5000, ...
        'MaxIterations', 1000);
    
    % Run optimization
    fprintf('Fitting extrema detection model...\n');
    [x_opt, neg_log_lik, exitflag] = fmincon(obj_fun, x0, [], [], [], [], lb, ub, [], options);
    
    % Package results
    params.bound = x_opt(1);
    params.t_nd_mean = x_opt(2);
    params.t_nd_std = x_opt(3);
    
    fprintf('\n=== Fitted Parameters ===\n');
    fprintf('Bound:         %.4f\n', params.bound);
    fprintf('T_ND mean:     %.4f s\n', params.t_nd_mean);
    fprintf('T_ND std:      %.4f s\n', params.t_nd_std);
    fprintf('Neg log-lik:   %.2f\n', neg_log_lik);
    fprintf('Exit flag:     %d\n', exitflag);
end


% Example usage
function example_usage()
    % Simulate some example data
    rng(42);
    
    % Single trial example
    n_timepoints = 1000;
    dt = 0.001;  % 1ms steps
    
    % Create a signal with drift + noise
    drift = 5.0;
    signal = normrnd(drift * sqrt(dt), sqrt(dt), n_timepoints, 1);
    
    % Simulate RT
    bound = 0.1;
    rt_observed = 0.35;  % 350ms
    
    % Calculate likelihood
    lik = extrema_detection_likelihood(signal, rt_observed, bound, dt);
    fprintf('Likelihood for single trial: %.6e\n', lik);
    
    % Multiple trials example
    n_trials = 50;
    signals = cell(n_trials, 1);
    for i = 1:n_trials
        signals{i} = normrnd(drift * sqrt(dt), sqrt(dt), n_timepoints, 1);
    end
    
    % Generate RTs from known parameters
    true_bound = 0.15;
    true_t_nd_mean = 0.32;
    true_t_nd_std = 0.06;
    
    rts = zeros(n_trials, 1);
    for i = 1:n_trials
        % Find first crossing
        crossing_idx = find(signals{i} > true_bound, 1);
        if isempty(crossing_idx)
            decision_time = n_timepoints * dt;
        else
            decision_time = (crossing_idx - 1) * dt;
        end
        % Add non-decision time
        rts(i) = decision_time + normrnd(true_t_nd_mean, true_t_nd_std);
    end
    
    fprintf('\n=== True Parameters ===\n');
    fprintf('Bound:         %.4f\n', true_bound);
    fprintf('T_ND mean:     %.4f s\n', true_t_nd_mean);
    fprintf('T_ND std:      %.4f s\n', true_t_nd_std);
    
    % Fit the model
    [params, neg_log_lik, exitflag] = fit_extrema_model(signals, rts, dt, ...
        'InitBound', 0.1, ...
        'InitTndMean', 0.3, ...
        'InitTndStd', 0.05);
end