function results = mv_reg(data, varargin)
    % Multivariate linear regression for freeze duration differences
    %
    %
    % INPUTS:
    %   data - table with variables:
    %       .diffRT              : RT difference (focal - template match) [frames]
    %       .diff_AcEv_before    : Diff in accumulated evidence before template
    %       .diff_AcEv_during    : Diff in accumulated evidence during template
    %       .condition           : Condition labels (will be treated as categorical)
    %
    % OPTIONS:
    %   'plots'       - logical, create diagnostic plots (default: true)
    %   'alpha'       - significance level (default: 0.05)
    %   'bootstrap'   - perform bootstrap CI estimation (default: false)
    %   'n_boot'      - number of bootstrap samples (default: 1000)
    %
    % OUTPUTS:
    %   results - comprehensive structure with all analysis results
    
    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'plots', true, @islogical);
    addParameter(p, 'alpha', 0.05, @isnumeric);
    addParameter(p, 'bootstrap', false, @islogical);
    addParameter(p, 'n_boot', 1000, @isnumeric);
    parse(p, varargin{:});
    opts = p.Results;
    
    
    % Data preprocessing
    fprintf('── Data Preprocessing ──\n');
    n_original = height(data);
    data = rmmissing(data);
    n_clean = height(data);
    fprintf('Original observations: %d\n', n_original);
    fprintf('After removing missing: %d\n', n_clean);
    if n_original > n_clean
        fprintf('⚠ Removed %d observations with missing data\n', n_original - n_clean);
    end
    fprintf('\n');
    
    % Descriptive statistics
    fprintf('── Descriptive Statistics ──\n');
    fprintf('diffRT:              M = %7.2f, SD = %6.2f, range = [%7.2f, %7.2f]\n', ...
        mean(data.diffRT), std(data.diffRT), min(data.diffRT), max(data.diffRT));
    fprintf('diff_AcEv_before:    M = %7.2f, SD = %6.2f, range = [%7.2f, %7.2f]\n', ...
        mean(data.diff_AcEv_before), std(data.diff_AcEv_before), ...
        min(data.diff_AcEv_before), max(data.diff_AcEv_before));
    fprintf('diff_AcEv_during:    M = %7.2f, SD = %6.2f, range = [%7.2f, %7.2f]\n', ...
        mean(data.diff_AcEv_during), std(data.diff_AcEv_during), ...
        min(data.diff_AcEv_during), max(data.diff_AcEv_during));
    
    % Check predictor correlation
    r_pred = corr(data.diff_AcEv_before, data.diff_AcEv_during);
    fprintf('\nCorrelation between predictors: r = %.3f', r_pred);
    if abs(r_pred) > 0.7
        fprintf(' ⚠ HIGH CORRELATION - multicollinearity concern!\n');
    else
        fprintf(' ✓\n');
    end
    fprintf('\n');
    
    % Fit the model
    fprintf('── Model Fitting ──\n');
    fprintf('Formula: diffRT ~ diff_AcEv_before + diff_AcEv_during \n\n');
    
    mdl = fitlm(data, 'diffRT ~ diff_AcEv_before + diff_AcEv_during');
    
    % Extract results
    results.model = mdl;
    results.data = data;
    results.n_obs = n_clean;
    results.coefficients = mdl.Coefficients;
    results.R2 = mdl.Rsquared.Ordinary;
    results.R2_adj = mdl.Rsquared.Adjusted;
    results.RMSE = mdl.RMSE;
    results.F_stat = mdl.anova.F(1);
    results.F_pval = mdl.anova.pValue(1);
    results.AIC = mdl.ModelCriterion.AIC;
    results.BIC = mdl.ModelCriterion.BIC;
    
    % Store key coefficients
    coef_names = mdl.CoefficientNames;
    idx_before = find(strcmp(coef_names, 'diff_AcEv_before'));
    idx_during = find(strcmp(coef_names, 'diff_AcEv_during'));
    
    results.beta_before = mdl.Coefficients.Estimate(idx_before);
    results.beta_during = mdl.Coefficients.Estimate(idx_during);
    results.se_before = mdl.Coefficients.SE(idx_before);
    results.se_during = mdl.Coefficients.SE(idx_during);
    results.p_before = mdl.Coefficients.pValue(idx_before);
    results.p_during = mdl.Coefficients.pValue(idx_during);
    results.t_before = mdl.Coefficients.tStat(idx_before);
    results.t_during = mdl.Coefficients.tStat(idx_during);
    
    % Confidence intervals
    CI = coefCI(mdl, 1-opts.alpha);
    results.CI_before = CI(idx_before, :);
    results.CI_during = CI(idx_during, :);
    
    % Display results
    disp(mdl);