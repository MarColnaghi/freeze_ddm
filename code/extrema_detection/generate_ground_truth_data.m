%% MULTIVARIATE REGRESSION ANALYSIS FOR RT AND ACCUMULATED EVIDENCE
% This package provides a complete statistical framework for analyzing how
% accumulated evidence before and during template presentation affects RT
% differences, while controlling for experimental conditions.
%
% Author: MATLAB Statistical Analysis Expert
% Date: October 2025

%% ========================================================================
% PART 1: GROUND TRUTH SIMULATION
% ========================================================================
% This simulates data with KNOWN effects so we can validate the analysis

function [data, ground_truth] = generate_ground_truth_data(n_obs, params)
    % Generate synthetic data with known ground truth parameters
    %
    % INPUTS:
    %   n_obs  - number of observations (default: 500)
    %   params - structure with ground truth parameters
    %
    % OUTPUTS:
    %   data         - table with simulated observations
    %   ground_truth - structure with true parameter values
    
    if nargin < 1, n_obs = 500; end
    if nargin < 2
        % Default ground truth parameters
        params.beta_intercept = 0;           % Baseline diffRT
        params.beta_before = -35;            % Effect of pre-template evidence (ms per unit)
        params.beta_during = -5;             % Effect of during-template evidence (small!)
        params.beta_condition2 = 15;         % Effect of condition 2 vs condition 1
        params.beta_condition3 = -10;        % Effect of condition 3 vs condition 1
        params.noise_sd = 25;                % Residual standard deviation
    end
    
    fprintf('=== GENERATING GROUND TRUTH DATA ===\n');
    fprintf('Sample size: %d observations\n', n_obs);
    fprintf('True parameters:\n');
    fprintf('  β₀ (intercept)           = %.2f ms\n', params.beta_intercept);
    fprintf('  β₁ (pre-template AcEv)   = %.2f ms/unit\n', params.beta_before);
    fprintf('  β₂ (during-template AcEv) = %.2f ms/unit\n', params.beta_during);
    fprintf('  β₃ (condition 2 effect)  = %.2f ms\n', params.beta_condition2);
    fprintf('  β₄ (condition 3 effect)  = %.2f ms\n', params.beta_condition3);
    fprintf('  Residual SD              = %.2f ms\n\n', params.noise_sd);
    
    % Generate predictor variables
    % Accumulated evidence differences (somewhat correlated, realistic)
    rho = 0.3; % correlation between before and during
    mu = [0, 0];
    sigma = [1, rho; rho, 1];
    AcEv = mvnrnd(mu, sigma, n_obs);
    
    diff_AcEv_before = AcEv(:, 1) * 2;      % Scale to reasonable range
    diff_AcEv_during = AcEv(:, 2) * 0.8;    % Smaller range (better matching)
    
    % Generate conditions (balanced design)
    condition = repmat([1; 2; 3], ceil(n_obs/3), 1);
    condition = condition(1:n_obs);
    condition = condition(randperm(n_obs));  % Shuffle
    
    % Create design matrix for conditions (dummy coding)
    X_cond2 = double(condition == 2);
    X_cond3 = double(condition == 3);
    
    % Generate response variable with TRUE effects
    diffRT = params.beta_intercept + ...
             params.beta_before * diff_AcEv_before + ...
             params.beta_during * diff_AcEv_during + ...
             params.beta_condition2 * X_cond2 + ...
             params.beta_condition3 * X_cond3 + ...
             randn(n_obs, 1) * params.noise_sd;
    
    % Assemble into table
    data = table(diffRT, diff_AcEv_before, diff_AcEv_during, condition, ...
                 'VariableNames', {'diffRT', 'diff_AcEv_before', ...
                                   'diff_AcEv_during', 'condition'});
    
    % Store ground truth
    ground_truth = params;
    ground_truth.correlation_before_during = rho;
    ground_truth.sample_size = n_obs;
    
    fprintf('Data generated successfully!\n');
    fprintf('Correlation between predictors: r = %.2f\n\n', rho);
end
