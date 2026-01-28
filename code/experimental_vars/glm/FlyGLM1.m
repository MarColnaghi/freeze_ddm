id_code = 'imm3_mob3_pc4';
paths_out = path_generator('folder', 'descriptive/fd_durs', 'bouts_id', id_code, 'imfirst', false);

% col = cmapper();
thresholds = define_thresholds;
thresholds.le_window_fl = [5 40];
thresholds.le_window_sl = [15 50];

bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);

paths_out = path_generator('folder', '/spontaneous_process', 'bouts_id', id_code, 'imfirst', false);
bouts_le = data_parser_new(bouts, 'period', 'loom', 'window', 'le', 'type', 'immobility', 'nloom', 2:20, 'min_dur', 30);

T = bouts_le;

Tmax = 10.5;
Ttrunc = 0.5;

% change flyID to Subject
VarNames = T.Properties.VariableNames;
VarNames{strcmp(VarNames,'fly')} = 'Subject';
T.Properties.VariableNames = VarNames;

Tmax = 10.5;
Ttrunc = 0.5;
% T.durations_s(T.durations_s>Tmax) = Tmax;
% truncate
T = T(T.durations_s>Ttrunc,:);

% apply log transform?
% T.durations_s = log(T.durations_s);

%% 2. PREPROCESSING: Z-SCORE EVERYTHING (Stores Mu and Sigma)
T_fit = T;
preds = {'nloom', 'avg_fs_1s', 'moving_flies', 'sloom', 'avg_sm','time_since_last', 'n_generated_freezes', 'cum_freeze_time'};
%time_since_last, n_generated_freezes, cum_freeze_time
scale_map = struct();

for p = 1:length(preds)
    % 1. Force Double
    raw_col = double(T.(preds{p}));
    z_col_name = [preds{p}];


    mu = mean(raw_col, 'omitnan');
    sigma = std(raw_col, 'omitnan');

    if sigma < 1e-12
        warning('Variable %s is constant! Z-score set to 0.', preds{p});
        z_vals = zeros(size(raw_col));
        sigma = 1;
    else
        z_vals = (raw_col - mu) / sigma;
    end
    % Store as vectors matching table height
    mu_vec = repmat(mu, height(T), 1);
    sigma_vec = repmat(sigma, height(T), 1);

    % Diagnostic Check
    z_vals(isnan(z_vals)) = 0;
    z_vals(isinf(z_vals)) = 0;

    T_fit.(z_col_name) = z_vals;

    % SAVE BOTH PARAMETERS
    scale_map.(preds{p}).Mu = mu_vec;
    scale_map.(preds{p}).Sigma = sigma_vec;
end

% --- C. GLMM ---
disp('Fitting GLMM...');
formula_mixed = ['durations_s ~ 1 + nloom + avg_fs_1s + moving_flies + sloom + avg_sm + time_since_last + n_generated_freezes + cum_freeze_time + avg_history_dur + (1 | Subject)'];% ...
   % ' (1 +  nloom + avg_fs_1s + moving_flies + sloom + avg_sm + time_since_last + n_generated_freezes| Subject)'];


glmm = fitglme(T_fit, formula_mixed, ...
    'Distribution', 'Gamma', ...
    'Link', 'log',...
    'CovariancePattern', 'Diagonal', ...
    'Exclude', []); % 'Exclude', 'none' forces it to fail loudly if NaNs persist

% normal GLM
disp('Fitting GLMM...');
formula_glm = 'durations_s ~ 1 + nloom + avg_fs_1s + moving_flies + sloom + avg_sm';


glm = fitglm(T_fit, formula_glm, ...
    'Distribution', 'Normal'); % 'Exclude', 'none' forces it to fail loudly if NaNs persist



% %% 4. EXTRACTING & UN-SCALING BETAS (RECOVERING "RAW" UNITS)
% 
% % Prepare Data containers for plotting
% params_to_plot = preds;
% titles_plot    = preds;
% % truths         = [beta_stim, beta_hist, beta_peia];
% 
% % Data Structures for plotting [SubjectID, Value]
% plot_data = struct(); 
% 
% for p = 1:length(params_to_plot)
%     p_name = params_to_plot{p};
% 
%     % --- 1. Get Sigmas for Un-scaling ---
%     if contains(p_name, ':') 
%         parts = split(p_name, ':');
%         term1 = strrep(parts{1}, '_Z', '');
%         term2 = strrep(parts{2}, '_Z', '');
%         % Access .Sigma specifically
%         sigma_vec = scale_map.(term1).Sigma .* scale_map.(term2).Sigma;
%     else
%         term = strrep(p_name, '_Z', '');
%         sigma_vec = scale_map.(term).Sigma;
%     end
% 
%     % Average Sigma per subject (to unscale subject estimates)
%     G_sub = findgroups(T_fit.Subject);
%     sigma_per_sub = splitapply(@mean, sigma_vec, G_sub);
%     % Average Sigma global (to unscale pooled/fixed estimates)
%     sigma_global = mean(sigma_vec);
% 
%     % % --- 2. Pooled Estimate ---
%     % b_pool = glm_pool.Coefficients.Estimate(p_name);
%     % ci_pool = glm_pool.coefCI; 
%     % ci_pool = ci_pool(strcmp(glm_pool.CoefficientNames, p_name), :);
%     % 
%     % plot_data(p).Pooled = b_pool / sigma_global;
%     % plot_data(p).PooledCI = ci_pool / sigma_global;
% 
% 
% 
%     % --- 4. GLMM BLUPs (Subject Specific) ---
%     % Fixed Effect
%     fixed_b = glmm.Coefficients.Estimate(find(strcmp(glmm.CoefficientNames, p_name)));
% 
%     % Random Effects (BLUPs)
%     [b_rand, re_names] = randomEffects(glmm);
% 
%     % Filter for Subject Level REs matching this parameter
%     % Note: fitglme names REs like 'Stim_Z'
%     is_target = strcmp(re_names.Name, p_name) & strcmp(re_names.Group, 'Subject');
% 
%     % Match BLUPs to Subjects
%     re_vals = b_rand(is_target);
%     re_subs = re_names.Level(is_target); 
% 
%     % Combine Fixed + Random, then Unscale
%     glmm_vals = nan(length(subs), 1);
% 
%     % Also get Population Estimate (Fixed Only)
%     plot_data(p).GLMM_Fixed = fixed_b / sigma_global;
% 
%     % We need to map the REs back to the sorted subject order 1..20
%     for i = 1:length(subs)
%         sub_str = string(subs(i));
% 
%         % Find RE for this sub (Compare string to string)
%         idx = string(re_subs) == sub_str;
% 
%         if any(idx)
%             total_beta_z = fixed_b + re_vals(idx);
%             glmm_vals(i) = total_beta_z / sigma_per_sub(i);
%         else
%             % If no RE (shrinkage to 0), just Fixed
%             glmm_vals(i) = fixed_b / sigma_per_sub(i);
%         end
%     end
%     plot_data(p).GLMM_BLUPs = glmm_vals;
% end
% 
% %% 5. VISUALIZATION: THE CATERPILLAR PLOT (Full Comparison)
% 
% figure('Color','w', 'Position', [50 50 1400 500]);
% 
% % Map plot indices to ParamsT column names
% param_map = {'Slope', 'BetaHist', 'BetaPEIA'};
% 
% for p = 1:3
%     subplot(1, 3, p); hold on;
%     d = plot_data(p);
%     x = 1:nSubjects;
% 
%     % --- RETRIEVE TRUE INDIVIDUAL VALUES ---
%     % Get the column name corresponding to this parameter
%     t_col = param_map{p};
%     % Calculate true mean per subject from the Generation Table (ParamsT)
%     % Ensure sorting matches the 1..20 index
%     true_sub_vals = splitapply(@mean, ParamsT.(t_col), findgroups(ParamsT.Subject));
%     true_empirical_mean = mean(true_sub_vals);
% 
%     % --- REFERENCE LINES ---
%     % 1. Theoretical Truth (Green Solid)
%     yline(truths(p), 'g-', 'LineWidth', 2, 'DisplayName', 'Theoretical Truth');
% 
%     % 2. Empirical Truth (Green Dashed) - NEW
%     yline(true_empirical_mean, 'g--', 'LineWidth', 1.5, ...
%         'DisplayName', 'True Empirical Mean');
% 
%     % 3. Pooled GLM (Blue Dashed)
%     yline(d.Pooled, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Pooled GLM');
% 
%     % 4. GLMM Fixed Effect (Red Dashed)
%     yline(d.GLMM_Fixed, 'r--', 'LineWidth', 1.5, 'DisplayName', 'GLMM Fixed Effect');
% 
%     % --- SCATTER PLOTS ---
%     % 5. Individual GLMs (Grey Triangles)
%     s1 = scatter(x, d.Indiv, 40, '^', 'MarkerEdgeColor', [0.5 0.5 0.5], ...
%         'MarkerFaceColor', [0.8 0.8 0.8], 'DisplayName', 'Individual GLM Fits');
% 
%     % 6. True Subject Betas (Green Diamonds) - NEW
%     s_true = scatter(x, true_sub_vals, 60, 'd', 'MarkerEdgeColor', [0 0.5 0], ...
%         'MarkerFaceColor', 'g', 'DisplayName', 'True Subject Beta');
% 
%     % 7. GLMM BLUPs (Black Dots)
%     s2 = scatter(x, d.GLMM_BLUPs, 40, 'k', 'filled', 'DisplayName', 'GLMM BLUPs');
% 
%     % 8. Connector Lines (True -> GLMM) - Visualizing Accuracy
%     % Connect the Green Diamond (Truth) to the Black Dot (Estimate)
%     plot([x; x], [true_sub_vals'; d.GLMM_BLUPs'], '-', 'Color', [0 0 0, 0.3], 'HandleVisibility', 'off');
% 
%     % --- OUTLIER PROTECTION (ROBUST ZOOM) ---
%     all_vals = [d.Indiv; d.GLMM_BLUPs; true_sub_vals];
%     center = median(all_vals, 'omitnan');
%     robust_sigma = 1.4826 * mad(all_vals, 1); 
%     k = 4; 
%     y_min = center - (k * robust_sigma);
%     y_max = center + (k * robust_sigma);
%     interesting_lines = [truths(p), d.GLMM_Fixed, d.Pooled, true_empirical_mean];
%     y_min = min([y_min, interesting_lines]);
%     y_max = max([y_max, interesting_lines]);
%     range_span = y_max - y_min;
%     if range_span == 0, range_span = 0.1; end 
%     ylim([y_min - 0.1*range_span, y_max + 0.1*range_span]);
%     % ---------------------------------------
% 
%     % Styling
%     xlabel('Subject ID'); ylabel('Beta (Un-scaled)');
%     title(titles_plot{p}); grid on; xlim([0 nSubjects+1]);
% 
%     if p == 2 
%         legend([s_true, s2, s1], 'Location', 'Best'); 
%     end
% 
%     % Text Stats
%     txt = sprintf('True Emp: %.2f\nGLMM Fix: %.2f', true_empirical_mean, d.GLMM_Fixed);
%     yl = ylim; 
%     text(1, yl(1) + 0.1*(yl(2)-yl(1)), txt, 'FontSize', 8, 'BackgroundColor', 'w');
% end
% sgtitle(['Recovery of Betas (Scaling: ' SCALING_METHOD ')']);



%% OLD PLOTTING CODE


% =========================================================================
% 1. DATA EXTRACTION AND PREPARATION
%    Extract all necessary data into a structured format first.
% =========================================================================
disp('Preparing data for plotting...');
model_results = struct(); % Initialize a struct to hold all results


subject_level_group_name = 'Subject';
fprintf('Plotting sigmas for subject-level group: ''%s''\n', subject_level_group_name);

full_model = glmm;

% Store basic results
model_results.betas = full_model.Coefficients.Estimate;
model_results.tStats = full_model.Coefficients.tStat;
model_results.betanames = full_model.CoefficientNames;
model_results.SE = full_model.Coefficients.SE;

% =========================================================================
% --- Dynamically extract sigmas for the subject-level random effects ---
% =========================================================================
model_results.sigma_map = containers.Map('KeyType', 'char', 'ValueType', 'double');

% 1. Get the list of all random effect grouping variables from the formula
re_group_names = full_model.Formula.GroupingVariableNames;

% 2. Find the index of our subject-level group.
%    cellfun applies a function to each cell of re_group_names. The function
%    (@(c) ismember(...)) looks *inside* each inner cell 'c' to see if it
%    contains our subject_level_group_name. This returns a logical array
%    (e.g., [0 1]), which `find` then uses to get the correct index.
re_index = find(cellfun(@(c) ismember(subject_level_group_name, c), re_group_names));

if ~isempty(re_index)
    % 3. Get the correct covariance matrix using the dynamically found index.
    cov_matrix = full_model.covarianceParameters{re_index};

    % 4. Calculate the standard deviations (sigmas) from this matrix.
    re_sigmas_raw = sqrt(diag(cov_matrix));

    % 5. Get the LinearFormula object for this specific random effect.
    re_formula_object = full_model.Formula.RELinearFormula{re_index};

    % 6. GET THE NAMES DIRECTLY FROM THE PROPERTIES YOU FOUND!
    %    The .TermNames property includes the intercept, which is what we need to
    %    match the full covariance matrix.
    re_coeff_names = re_formula_object.TermNames;

    % 7. Store the sigmas in a map. The order of sigmas from the diagonal
    %    of the covariance matrix directly corresponds to the order of TermNames.
    for j = 1:length(re_coeff_names)
        model_results.sigma_map(re_coeff_names{j}) = re_sigmas_raw(j);
    end
end

% --- Find the common predictors and groups to plot ---
% (The rest of the plotting code from here down is correct and does not need to change)
common_betanames = model_results.betanames;
% --- uncomment to remove the intercept from the plots
% common_betanames = common_betanames(~strcmp(common_betanames, '(Intercept)'));

common_group_names = intersect(preds, preds, 'stable'); % Assumes var_group_names is the same

n_betas_toplot = length(common_betanames);
n_groups_toplot = length(common_group_names);

% =========================================================================
% 2. PLOTTING USING TILEDLAYOUT (Horizontal orientation)
% =========================================================================
disp('Plotting GLMM results with tiledlayout (horizontal)...');

ff = figure('Color','w','Name','GLMM results');
ff.Units = 'centimeters';
ff.Position = [2 2.5 16 8];

% Create a 1x2 tiled layout
t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% Common y-axis positions and labels (horizontal bars)
plot_yticks = 1:n_betas_toplot;
plot_labels = cellfun(@(x) strrep(x, '_', '\_'), common_betanames, 'UniformOutput', false);
plot_labels = cellfun(@(x) strrep(x, '\_\_', '*'), plot_labels, 'UniformOutput', false);

col = get(0,'defaultAxesColorOrder');

% --- PLOT 1: Betas (horizontal) ---
nexttile(1); hold on;

% Reference vertical line at 0
xline(0, 'k--');

for i = 1:n_betas_toplot
    pred_name = common_betanames{i};
    
    beta_idx = find(strcmp(model_results.betanames, pred_name));
    beta_val = model_results.betas(beta_idx);
    se_val = model_results.SE(beta_idx);

    if isKey(model_results.sigma_map, pred_name)
        sigma_val = model_results.sigma_map(pred_name);
    else
        sigma_val = 0;
    end

    total_sigma = sqrt(se_val.^2 + sigma_val.^2);

    % Horizontal point
    plot(beta_val, i, '.', 'MarkerSize', 17, 'Color', col(1,:));

    % Error bars (horizontal)
    if total_sigma > 0
        errorbar(beta_val, i, total_sigma, 'horizontal', 'Color', col(1,:));
        % Optional dotted line for random-effect sigma only
        errorbar(beta_val, i+0.05, sigma_val, 'horizontal', 'Color', col(1,:), 'LineStyle', ':');
    end
end

ylim([0, n_betas_toplot+1]);
yticks(plot_yticks); yticklabels(plot_labels);
xlabel('Betas'); grid on;

% --- PLOT 2: t-values (horizontal) ---
nexttile(2); hold on;

% Reference vertical lines
xline(0, 'k--');   % zero line
xline(2, 'k:');    % t = +2
xline(-2, 'k:');   % t = -2

for i = 1:n_betas_toplot
    pred_name = common_betanames{i};
    tStat_idx = find(strcmp(model_results.betanames, pred_name));
    tStat_val = model_results.tStats(tStat_idx);

    % Horizontal point
    plot(tStat_val, i, '.', 'MarkerSize', 17, 'Color', col(1,:));
end

ylim([0, n_betas_toplot+1]);
yticks(plot_yticks); yticklabels(plot_labels);
xlabel('t-values'); grid on;

% Optional: titles
title(t.Children(2), 'Betas');    % first tile
title(t.Children(1), 't-values'); % second tile

