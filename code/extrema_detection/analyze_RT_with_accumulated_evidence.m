function results = analyze_RT_with_accumulated_evidence(data, varargin)
    % Multivariate linear regression for RT differences
    %
    % USAGE:
    %   results = analyze_RT_with_accumulated_evidence(data)
    %   results = analyze_RT_with_accumulated_evidence(data, 'plots', true)
    %
    % INPUTS:
    %   data - table with variables:
    %       .diffRT              : RT difference (focal - template match) [ms]
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
    
    fprintf('╔════════════════════════════════════════════════════════════╗\n');
    fprintf('║   MULTIVARIATE REGRESSION ANALYSIS: RT vs AcEv            ║\n');
    fprintf('╚════════════════════════════════════════════════════════════╝\n\n');
    
    % Convert to table if needed
    if isstruct(data)
        data = struct2table(data);
    end
    
    % Ensure condition is categorical
    if isnumeric(data.condition)
        data.condition = categorical(data.condition);
    end
    
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
    fprintf('Formula: diffRT ~ diff_AcEv_before + diff_AcEv_during + condition\n\n');
    
    mdl = fitlm(data, 'diffRT ~ diff_AcEv_before + diff_AcEv_during + condition');
    
    % Extract results
    results.model = mdl;
    results.data = data;
    results.n_obs = n_clean;
    results.coefficients = mdl.Coefficients;
    results.R2 = mdl.Rsquared.Ordinary;
    results.R2_adj = mdl.Rsquared.Adjusted;
    results.RMSE = mdl.RMSE;
%     results.F_stat = mdl.anova.F(1);
%     results.F_pval = mdl.anova.pValue(1);
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
    
    % Bootstrap confidence intervals (if requested)
    if opts.bootstrap
        fprintf('\n── Bootstrap Analysis ──\n');
        fprintf('Computing %d bootstrap samples...\n', opts.n_boot);
        boot_results = bootstrp(opts.n_boot, @(idx) fit_and_extract(data(idx,:)), ...
                                (1:height(data))');
        results.bootstrap.beta_before = boot_results(:,1);
        results.bootstrap.beta_during = boot_results(:,2);
        results.bootstrap.CI_before = prctile(boot_results(:,1), [opts.alpha/2, 1-opts.alpha/2]*100);
        results.bootstrap.CI_during = prctile(boot_results(:,2), [opts.alpha/2, 1-opts.alpha/2]*100);
        fprintf('Bootstrap complete!\n\n');
    end
    
    % Model diagnostics
    fprintf('── Model Diagnostics ──\n');
    fprintf('R² = %.4f (%.1f%% variance explained)\n', results.R2, results.R2*100);
    fprintf('Adjusted R² = %.4f\n', results.R2_adj);
    fprintf('RMSE = %.2f ms\n', results.RMSE);
%     fprintf('F(%d, %d) = %.2f, p %s\n', mdl.DFE, mdl.NumObservations-mdl.NumCoefficients, ...
%             results.F_stat, format_pvalue(results.F_pval));
    fprintf('AIC = %.1f, BIC = %.1f\n\n', results.AIC, results.BIC);
    
    % Check assumptions
    fprintf('── Assumption Checks ──\n');
    [h_sw, p_sw] = swtest(mdl.Residuals.Raw);
    fprintf('Shapiro-Wilk normality test: W = %.4f, p %s', 1-h_sw, format_pvalue(p_sw));
    if h_sw
        fprintf(' ⚠ Residuals non-normal\n');
    else
        fprintf(' ✓\n');
    end
    
    % Durbin-Watson for autocorrelation
    dw = dwtest(mdl);
    fprintf('Durbin-Watson statistic: DW = %.3f', dw);
    if dw < 1.5 || dw > 2.5
        fprintf(' ⚠ Potential autocorrelation\n');
    else
        fprintf(' ✓\n');
    end
    
    % VIF for multicollinearity
    if exist('vif', 'file')
        vif_vals = diag(inv(corrcoef(mdl.Variables{:, mdl.PredictorNames})));
        max_vif = max(vif_vals);
        fprintf('Max VIF: %.2f', max_vif);
        if max_vif > 5
            fprintf(' ⚠ High multicollinearity\n');
        else
            fprintf(' ✓\n');
        end
    end
    fprintf('\n');
    
    % Key findings
    print_key_findings(results, opts.alpha);
    
    % Create plots
    if opts.plots
        create_comprehensive_plots(results, opts.alpha);
    end
    
    % Final interpretation
    print_interpretation(results, opts.alpha);
end

%% ========================================================================
% PART 3: HELPER FUNCTIONS
% ========================================================================

function betas = fit_and_extract(data)
    % Helper for bootstrap: fit model and extract key coefficients
    mdl = fitlm(data, 'diffRT ~ diff_AcEv_before + diff_AcEv_during + condition', ...
                'Verbose', false);
    coef_names = mdl.CoefficientNames;
    idx_before = find(strcmp(coef_names, 'diff_AcEv_before'));
    idx_during = find(strcmp(coef_names, 'diff_AcEv_during'));
    betas = [mdl.Coefficients.Estimate(idx_before), ...
             mdl.Coefficients.Estimate(idx_during)];
end

function str = format_pvalue(p)
    % Format p-value with appropriate precision
    if p < 0.001
        str = '< 0.001 ***';
    elseif p < 0.01
        str = sprintf('= %.3f **', p);
    elseif p < 0.05
        str = sprintf('= %.3f *', p);
    else
        str = sprintf('= %.3f n.s.', p);
    end
end

function stars = get_stars(p)
    % Get significance stars
    if p < 0.001
        stars = '***';
    elseif p < 0.01
        stars = '**';
    elseif p < 0.05
        stars = '*';
    else
        stars = 'n.s.';
    end
end

function print_key_findings(results, alpha)
    fprintf('╔════════════════════════════════════════════════════════════╗\n');
    fprintf('║                    KEY FINDINGS                            ║\n');
    fprintf('╚════════════════════════════════════════════════════════════╝\n\n');
    
    fprintf('┌─ PREDICTOR 1: Pre-Template Accumulated Evidence ─────────┐\n');
    fprintf('│  β = %7.2f ms/unit  (SE = %.2f)                         │\n', ...
            results.beta_before, results.se_before);
    fprintf('│  95%% CI: [%7.2f, %7.2f]                                 │\n', ...
            results.CI_before(1), results.CI_before(2));
    fprintf('│  t(%d) = %6.2f, p %s                           │\n', ...
            results.model.DFE, results.t_before, format_pvalue(results.p_before));
    fprintf('└───────────────────────────────────────────────────────────┘\n\n');
    
    fprintf('┌─ PREDICTOR 2: During-Template Accumulated Evidence ──────┐\n');
    fprintf('│  β = %7.2f ms/unit  (SE = %.2f)                         │\n', ...
            results.beta_during, results.se_during);
    fprintf('│  95%% CI: [%7.2f, %7.2f]                                 │\n', ...
            results.CI_during(1), results.CI_during(2));
    fprintf('│  t(%d) = %6.2f, p %s                           │\n', ...
            results.model.DFE, results.t_during, format_pvalue(results.p_during));
    fprintf('└───────────────────────────────────────────────────────────┘\n\n');
    
    % Effect size comparison
    ratio = abs(results.beta_before) / abs(results.beta_during);
    fprintf('Effect size ratio: |β₁|/|β₂| = %.2f\n', ratio);
    fprintf('(Pre-template effect is %.1fx %s than during-template effect)\n\n', ...
            ratio, iif(ratio > 1, 'LARGER', 'smaller'));
end

function print_interpretation(results, alpha)
    fprintf('╔════════════════════════════════════════════════════════════╗\n');
    fprintf('║                   INTERPRETATION                           ║\n');
    fprintf('╚════════════════════════════════════════════════════════════╝\n\n');
    
    % Check for ideal pattern
    pre_sig = results.p_before < alpha;
    pre_neg = results.beta_before < 0;
    during_small = abs(results.beta_during) < abs(results.beta_before) / 2;
    during_nonsig = results.p_during >= alpha;
    
    if pre_sig && pre_neg
        fprintf('✓ POSITIVE FINDING: Pre-template accumulated evidence has a\n');
        fprintf('  significant negative effect on RT (p %s).\n', format_pvalue(results.p_before));
        fprintf('  → More evidence BEFORE template → FASTER responses\n\n');
    else
        fprintf('✗ Pre-template effect is not significant (p %s)\n\n', ...
                format_pvalue(results.p_before));
    end
    
    if during_nonsig || during_small
        fprintf('✓ POSITIVE FINDING: During-template evidence effect is %s\n', ...
                iif(during_nonsig, 'non-significant', 'small'));
        fprintf('  (β = %.2f, p %s)\n', results.beta_during, ...
                format_pvalue(results.p_during));
        fprintf('  → Template matching appears adequate\n');
        fprintf('  → Pre-template effect is NOT explained away by template differences\n\n');
    else
        fprintf('⚠ WARNING: During-template evidence has substantial effect\n');
        fprintf('  (β = %.2f, p %s)\n', results.beta_during, ...
                format_pvalue(results.p_during));
        fprintf('  → Template matching may not be perfect\n');
        fprintf('  → Part of the effect may be due to template mismatches\n\n');
    end
    
    % Overall conclusion
    fprintf('── OVERALL CONCLUSION ──\n');
    if pre_sig && pre_neg && (during_nonsig || during_small)
        fprintf('✓✓✓ STRONG EVIDENCE for history-dependent RT modulation!\n');
        fprintf('    Accumulated evidence before template onset has a unique,\n');
        fprintf('    robust effect on reaction times that cannot be explained\n');
        fprintf('    by template matching artifacts or condition differences.\n');
    elseif pre_sig && pre_neg
        fprintf('✓ MODERATE EVIDENCE for history-dependent RT modulation.\n');
        fprintf('  Pre-template effect exists but during-template evidence also\n');
        fprintf('  contributes. Consider improving template matching.\n');
    else
        fprintf('✗ WEAK/NO EVIDENCE for history-dependent effect.\n');
        fprintf('  Pre-template accumulated evidence does not significantly\n');
        fprintf('  predict RT differences after controlling for other factors.\n');
    end
    fprintf('\n');
end

function result = iif(condition, true_val, false_val)
    % Inline if function
    if condition
        result = true_val;
    else
        result = false_val;
    end
end

function create_comprehensive_plots(results, alpha)
    % Create comprehensive diagnostic and results plots
    
    figure('Position', [50 50 1400 900], 'Color', 'w');
    
    % Extract data
    mdl = results.model;
    data = results.data;
    fitted = mdl.Fitted;
    resid = mdl.Residuals.Raw;
    
    % 1. Main effects plot - Pre-template
    subplot(3, 4, 1);
    scatter(data.diff_AcEv_before, data.diffRT, 30, [0.3 0.3 0.3], ...
            'filled', 'MarkerFaceAlpha', 0.3);
    hold on;
    % Add regression line (marginal)
    x_range = [min(data.diff_AcEv_before), max(data.diff_AcEv_before)];
    y_pred = results.beta_before * x_range + mean(data.diffRT - results.beta_before * data.diff_AcEv_before);
    plot(x_range, y_pred, 'r-', 'LineWidth', 2.5);
    xlabel('Diff AcEv BEFORE template', 'FontSize', 10);
    ylabel('diffRT (ms)', 'FontSize', 10);
    title(sprintf('Pre-Template Effect\nβ = %.2f %s', results.beta_before, ...
                  get_stars(results.p_before)), 'FontSize', 11);
    grid on; box on;
    xlim([-200 200])
    % 2. Main effects plot - During-template
    subplot(3, 4, 2);
    scatter(data.diff_AcEv_during, data.diffRT, 30, [0.3 0.3 0.3], ...
            'filled', 'MarkerFaceAlpha', 0.3);
    hold on;
    x_range = [min(data.diff_AcEv_during), max(data.diff_AcEv_during)];
    y_pred = results.beta_during * x_range + mean(data.diffRT - results.beta_during * data.diff_AcEv_during);
    plot(x_range, y_pred, 'b-', 'LineWidth', 2.5);
    xlabel('Diff AcEv DURING template', 'FontSize', 10);
    ylabel('diffRT (ms)', 'FontSize', 10);
    title(sprintf('During-Template Effect\nβ = %.2f %s', results.beta_during, ...
                  get_stars(results.p_during)), 'FontSize', 11);
    grid on; box on;
    xlim([-200 200])

    % 3. Predicted vs Actual
    subplot(3, 4, 3);
    scatter(data.diffRT, fitted, 30, [0.2 0.4 0.7], 'filled', 'MarkerFaceAlpha', 0.4);
    hold on;
    lims = [min([data.diffRT; fitted]), max([data.diffRT; fitted])];
    plot(lims, lims, 'k--', 'LineWidth', 1.5);
    xlabel('Actual diffRT (ms)', 'FontSize', 10);
    ylabel('Predicted diffRT (ms)', 'FontSize', 10);
    title(sprintf('Model Fit\nR² = %.3f', results.R2), 'FontSize', 11);
    grid on; box on;
    axis equal;
    xlim(lims); ylim(lims);
    
    % 4. Coefficient plot with CIs
    subplot(3, 4, 4);
    coef_idx = [find(strcmp(mdl.CoefficientNames, 'diff_AcEv_before')), ...
                find(strcmp(mdl.CoefficientNames, 'diff_AcEv_during'))];
    betas = [results.beta_before; results.beta_during];
    CIs = [results.CI_before; results.CI_during];
    errorbar(1:2, betas, betas - CIs(:,1), CIs(:,2) - betas, ...
             'o', 'MarkerSize', 10, 'MarkerFaceColor', [0.2 0.5 0.8], ...
             'LineWidth', 2, 'Color', [0.2 0.5 0.8], 'CapSize', 15);
    hold on;
    yline(0, 'k--', 'LineWidth', 1.5);
    xlim([0.5 2.5]);
    xticks(1:2);
    xticklabels({'Pre-Template', 'During-Template'});
    ylabel('Coefficient (ms/unit)', 'FontSize', 10);
    title('Effect Sizes with 95% CI', 'FontSize', 11);
    grid on; box on;
    
    % 5. Residuals vs Fitted
    subplot(3, 4, 5);
    scatter(fitted, resid, 30, [0.6 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.4);
    hold on;
    yline(0, 'k--', 'LineWidth', 1.5);
    % Add smoothed trend
    [f_sorted, idx] = sort(fitted);
    r_sorted = resid(idx);
    window = round(length(f_sorted) / 10);
    r_smooth = movmean(r_sorted, window);
    plot(f_sorted, r_smooth, 'b-', 'LineWidth', 2);
    xlabel('Fitted values', 'FontSize', 10);
    ylabel('Residuals', 'FontSize', 10);
    title('Residuals vs Fitted', 'FontSize', 11);
    grid on; box on;
    
    % 6. Q-Q Plot
    subplot(3, 4, 6);
    qqplot(resid);
    title('Normal Q-Q Plot', 'FontSize', 11);
    xlabel('Theoretical Quantiles', 'FontSize', 10);
    ylabel('Sample Quantiles', 'FontSize', 10);
    grid on;
    
    % 7. Histogram of residuals
    subplot(3, 4, 7);
    histogram(resid, 30, 'Normalization', 'pdf', 'FaceColor', [0.7 0.7 0.7], ...
              'EdgeColor', 'none', 'FaceAlpha', 0.7);
    hold on;
    x_norm = linspace(min(resid), max(resid), 100);
    y_norm = normpdf(x_norm, mean(resid), std(resid));
    plot(x_norm, y_norm, 'r-', 'LineWidth', 2);
    xlabel('Residuals', 'FontSize', 10);
    ylabel('Density', 'FontSize', 10);
    title('Residual Distribution', 'FontSize', 11);
    grid on; box on;
    legend('Observed', 'Normal', 'Location', 'best');
    
    % 8. Scale-Location plot
    subplot(3, 4, 8);
    scatter(fitted, sqrt(abs(resid)), 30, [0.4 0.6 0.4], ...
            'filled', 'MarkerFaceAlpha', 0.4);
    hold on;
    [f_sorted, idx] = sort(fitted);
    r_sorted = sqrt(abs(resid(idx)));
    r_smooth = movmean(r_sorted, window);
    plot(f_sorted, r_smooth, 'r-', 'LineWidth', 2);
    xlabel('Fitted values', 'FontSize', 10);
    ylabel('√|Residuals|', 'FontSize', 10);
    title('Scale-Location Plot', 'FontSize', 11);
    grid on; box on;
    
    % 9. Partial residual plot - Pre-template
    subplot(3, 4, 9);
    % Compute partial residuals for pre-template predictor
    idx_before = find(strcmp(mdl.CoefficientNames, 'diff_AcEv_before'));
    partial_resid = resid + results.beta_before * data.diff_AcEv_before;
    scatter(data.diff_AcEv_before, partial_resid, 30, [0.8 0.3 0.3], ...
            'filled', 'MarkerFaceAlpha', 0.4);
    hold on;
    x_range = [min(data.diff_AcEv_before), max(data.diff_AcEv_before)];
    y_pred = results.beta_before * x_range;
    plot(x_range, y_pred, 'k-', 'LineWidth', 2.5);
    xlabel('Diff AcEv DURING', 'FontSize', 10);
    ylabel('Partial Residuals', 'FontSize', 10);
    title('Partial: During-Template', 'FontSize', 11);
    grid on; box on;
    
    % 11. Cook's Distance
    subplot(3, 4, 11);
    cooks_d = mdl.Diagnostics.CooksDistance;
    stem(cooks_d, 'Marker', 'o', 'MarkerSize', 4, 'Color', [0.5 0.2 0.5], ...
         'MarkerFaceColor', [0.7 0.4 0.7]);
    hold on;
    thresh = 4 / results.n_obs;
    yline(thresh, 'r--', 'LineWidth', 1.5);
    xlabel('Observation', 'FontSize', 10);
    ylabel("Cook's Distance", 'FontSize', 10);
    title('Influential Points', 'FontSize', 11);
    grid on; box on;
    
    % 12. Leverage vs Residuals
    subplot(3, 4, 12);
    leverage = mdl.Diagnostics.Leverage;
    scatter(leverage, resid, 30, cooks_d, 'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
    yline(0, 'k--', 'LineWidth', 1.5);
    colorbar;
    xlabel('Leverage', 'FontSize', 10);
    ylabel('Residuals', 'FontSize', 10);
    title("Leverage vs Residuals (color = Cook's D)", 'FontSize', 11);
    grid on; box on;
    
    sgtitle('Comprehensive Regression Diagnostics', 'FontSize', 14, 'FontWeight', 'bold');
end

%% ========================================================================
% PART 4: VALIDATION FUNCTION - TEST RECOVERY OF GROUND TRUTH
% ========================================================================

function validation_results = validate_analysis_with_ground_truth()
    % Run complete validation: generate data with known parameters,
    % analyze it, and check if we recover the true parameters
    
    fprintf('\n');
    fprintf('════════════════════════════════════════════════════════════\n');
    fprintf('       GROUND TRUTH VALIDATION EXPERIMENT                   \n');
    fprintf('════════════════════════════════════════════════════════════\n\n');
    
    % Set ground truth parameters
    params.beta_intercept = 5;
    params.beta_before = -35;        % Strong negative effect
    params.beta_during = -5;         % Weak negative effect  
    params.beta_condition2 = 15;
    params.beta_condition3 = -10;
    params.noise_sd = 25;
    
    % Generate data
    [data, ground_truth] = generate_ground_truth_data(500, params);
    
    % Run analysis
    results = analyze_RT_with_accumulated_evidence(data, 'plots', true);
    
    % Compare recovered vs true parameters
    fprintf('\n');
    fprintf('════════════════════════════════════════════════════════════\n');
    fprintf('       PARAMETER RECOVERY ASSESSMENT                        \n');
    fprintf('════════════════════════════════════════════════════════════\n\n');
    
    fprintf('┌─────────────────────────────────────────────────────────┐\n');
    fprintf('│  Parameter          True Value   Recovered    Error     │\n');
    fprintf('├─────────────────────────────────────────────────────────┤\n');
    
    % Pre-template effect
    error_before = results.beta_before - ground_truth.beta_before;
    pct_error_before = abs(error_before / ground_truth.beta_before * 100);
    fprintf('│  β (pre-template)   %8.2f     %8.2f    %7.2f    │\n', ...
            ground_truth.beta_before, results.beta_before, error_before);
    
    % During-template effect
    error_during = results.beta_during - ground_truth.beta_during;
    pct_error_during = abs(error_during / ground_truth.beta_during * 100);
    fprintf('│  β (during-templ)   %8.2f     %8.2f    %7.2f    │\n', ...
            ground_truth.beta_during, results.beta_during, error_during);
    
    fprintf('└─────────────────────────────────────────────────────────┘\n\n');
    
    % Check if CIs contain true values
    ci_contains_before = ground_truth.beta_before >= results.CI_before(1) && ...
                         ground_truth.beta_before <= results.CI_before(2);
    ci_contains_during = ground_truth.beta_during >= results.CI_during(1) && ...
                         ground_truth.beta_during <= results.CI_during(2);
    
    fprintf('── Confidence Interval Coverage ──\n');
    fprintf('Pre-template:  95%% CI [%.2f, %.2f]  %s true value (%.2f)\n', ...
            results.CI_before(1), results.CI_before(2), ...
            iif(ci_contains_before, 'CONTAINS ✓', 'MISSES ✗'), ...
            ground_truth.beta_before);
    fprintf('During-template: 95%% CI [%.2f, %.2f]  %s true value (%.2f)\n\n', ...
            results.CI_during(1), results.CI_during(2), ...
            iif(ci_contains_during, 'CONTAINS ✓', 'MISSES ✗'), ...
            ground_truth.beta_during);
    
    % Statistical power check
    fprintf('── Statistical Power ──\n');
    fprintf('Pre-template effect:   p %s  %s\n', ...
            format_pvalue(results.p_before), ...
            iif(results.p_before < 0.05, '(Correctly detected ✓)', '(Missed! ✗)'));
    fprintf('During-template effect: p %s  %s\n\n', ...
            format_pvalue(results.p_during), ...
            iif(results.p_during < 0.05, '(Correctly detected ✓)', '(Missed! ✗)'));
    
    % Overall validation assessment
    fprintf('── VALIDATION SUMMARY ──\n');
    recovery_good = pct_error_before < 20 && pct_error_during < 30;
    ci_good = ci_contains_before && ci_contains_during;
    power_good = results.p_before < 0.05 && results.p_during < 0.05;
    
    if recovery_good && ci_good && power_good
        fprintf('✓✓✓ EXCELLENT: Analysis successfully recovers ground truth!\n');
        fprintf('    - Parameter estimates are accurate\n');
        fprintf('    - Confidence intervals have proper coverage\n');
        fprintf('    - Statistical power is adequate\n');
    elseif recovery_good && ci_good
        fprintf('✓✓ GOOD: Analysis recovers parameters accurately\n');
        fprintf('   (Some power issues with weaker effects)\n');
    elseif recovery_good
        fprintf('✓ ACCEPTABLE: Parameters reasonably well estimated\n');
        fprintf('  (Some CI coverage or power issues)\n');
    else
        fprintf('✗ ISSUES DETECTED in parameter recovery\n');
        fprintf('  Review model specification or increase sample size\n');
    end
    fprintf('\n');
    
    % Store validation results
    validation_results.ground_truth = ground_truth;
    validation_results.recovered = results;
    validation_results.error_before = error_before;
    validation_results.error_during = error_during;
    validation_results.ci_coverage = [ci_contains_before, ci_contains_during];
    validation_results.power_adequate = [results.p_before < 0.05, results.p_during < 0.05];
end

%% ========================================================================
% PART 5: POWER ANALYSIS
% ========================================================================

function power_results = perform_power_analysis(effect_sizes, sample_sizes, n_sim)
    % Monte Carlo power analysis for different effect sizes and sample sizes
    %
    % INPUTS:
    %   effect_sizes  - vector of beta values to test (e.g., -20:5:-40)
    %   sample_sizes  - vector of sample sizes (e.g., 100:100:500)
    %   n_sim         - number of simulations per combination (default: 1000)
    
    if nargin < 3, n_sim = 1000; end
    
    fprintf('\n');
    fprintf('════════════════════════════════════════════════════════════\n');
    fprintf('              POWER ANALYSIS (Monte Carlo)                  \n');
    fprintf('════════════════════════════════════════════════════════════\n\n');
    fprintf('Running %d simulations per condition...\n', n_sim);
    fprintf('This may take a few minutes...\n\n');
    
    n_effects = length(effect_sizes);
    n_samples = length(sample_sizes);
    power_matrix = zeros(n_effects, n_samples);
    
    % Progress tracking
    total_conditions = n_effects * n_samples;
    condition_count = 0;
    
    for i = 1:n_effects
        for j = 1:n_samples
            % Set parameters
            params.beta_intercept = 0;
            params.beta_before = effect_sizes(i);
            params.beta_during = -5;
            params.beta_condition2 = 10;
            params.beta_condition3 = -10;
            params.noise_sd = 25;
            
            % Run simulations
            p_values = zeros(n_sim, 1);
            for sim = 1:n_sim
                [data, ~] = generate_ground_truth_data(sample_sizes(j), params);
                mdl = fitlm(data, 'diffRT ~ diff_AcEv_before + diff_AcEv_during + condition', ...
                           'Verbose', false);
                idx = find(strcmp(mdl.CoefficientNames, 'diff_AcEv_before'));
                p_values(sim) = mdl.Coefficients.pValue(idx);
            end
            
            % Calculate power (proportion of significant results)
            power_matrix(i, j) = mean(p_values < 0.05);
            
            % Update progress
            condition_count = condition_count + 1;
            if mod(condition_count, 5) == 0
                fprintf('Progress: %d/%d conditions complete (%.1f%%)\n', ...
                        condition_count, total_conditions, ...
                        condition_count/total_conditions*100);
            end
        end
    end
    
    fprintf('\nPower analysis complete!\n\n');
    
    % Create power curve plot
    figure('Position', [100 100 1000 600], 'Color', 'w');
    
    subplot(1, 2, 1);
    imagesc(sample_sizes, effect_sizes, power_matrix);
    colorbar;
    colormap(hot);
    caxis([0 1]);
    xlabel('Sample Size', 'FontSize', 12);
    ylabel('Effect Size (β)', 'FontSize', 12);
    title('Statistical Power', 'FontSize', 13, 'FontWeight', 'bold');
    hold on;
    % Add 0.80 power contour
    contour(sample_sizes, effect_sizes, power_matrix, [0.8 0.8], 'b-', 'LineWidth', 2);
    set(gca, 'YDir', 'normal');
    
    subplot(1, 2, 2);
    colors = lines(n_effects);
    for i = 1:n_effects
        plot(sample_sizes, power_matrix(i, :), 'o-', 'LineWidth', 2, ...
             'MarkerSize', 8, 'Color', colors(i, :), ...
             'DisplayName', sprintf('β = %.0f', effect_sizes(i)));
        hold on;
    end
    yline(0.8, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Target power (0.80)');
    xlabel('Sample Size', 'FontSize', 12);
    ylabel('Statistical Power', 'FontSize', 12);
    title('Power Curves', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'southeast');
    grid on;
    ylim([0 1]);
    
    sgtitle('Power Analysis: Pre-Template Evidence Effect', ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    % Store results
    power_results.effect_sizes = effect_sizes;
    power_results.sample_sizes = sample_sizes;
    power_results.power_matrix = power_matrix;
    power_results.n_simulations = n_sim;
    
    % Print recommendations
    fprintf('── SAMPLE SIZE RECOMMENDATIONS ──\n');
    for i = 1:n_effects
        idx_80 = find(power_matrix(i, :) >= 0.80, 1, 'first');
        if ~isempty(idx_80)
            fprintf('For β = %6.1f: Need n ≥ %d for 80%% power\n', ...
                    effect_sizes(i), sample_sizes(idx_80));
        else
            fprintf('For β = %6.1f: Need n > %d for 80%% power\n', ...
                    effect_sizes(i), sample_sizes(end));
        end
    end
    fprintf('\n');
end
