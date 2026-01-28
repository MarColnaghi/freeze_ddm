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

% Apply truncation
% T = T(T.durations_s > Ttrunc, :);
% T.durations_s(T.durations_s > Tmax) = Tmax;

pct_censored = sum(T.durations_s > Tmax) / height(T) * 100;
disp(['Censored Data: ' num2str(pct_censored) '%']);

% apply log transform?
% T.durations_s = log(T.durations_s);

%% 2. PREPROCESSING: Z-SCORE EVERYTHING
T_fit = T;

% 1. Fix Response Scale (Use Milliseconds to prevent "Badly Scaled" warning)
T_fit.durations_ms = T_fit.durations_s * 1000;

% 2. Define ALL predictors (Added 'avg_ss' here)
preds = {'nloom', 'avg_fs_1s', 'moving_flies', 'sloom', 'avg_sm', 'avg_ss', ...
         'time_since_last', 'n_generated_freezes', 'cum_freeze_time'};

scale_map = struct();

for p = 1:length(preds)
    col_name = preds{p};
    
    % Force Double
    raw_col = double(T.(col_name));
    
    mu = mean(raw_col, 'omitnan');
    sigma = std(raw_col, 'omitnan');
    
    % Handle constants
    if sigma < 1e-12
        z_vals = zeros(size(raw_col));
        sigma = 1; 
    else
        z_vals = (raw_col - mu) / sigma;
    end
    
    % OVERWRITE column in T_fit with Z-scored data
    T_fit.(col_name) = z_vals;
    
    % Store parameters for later
    scale_map.(col_name).Mu = mu;       % Save scalar (efficient)
    scale_map.(col_name).Sigma = sigma; % Save scalar
end

%% 3. FITTING GLMM

disp('Fitting GLMM...');

% UPDATED FORMULA: 
% 1. Uses 'durations_ms' (stable scaling)
% 2. Includes Random Effects '(1|fly)' (required for fitglme)
formula_mixed = 'durations_s ~ 1 + nloom + sloom + avg_sm + avg_fs_1s + moving_flies + cum_freeze_time + (1 | fly)';

glmm = fitglme(T_fit, formula_mixed, ...
    'Distribution', 'InverseGaussian', ...
    'Link', 'log', ...       
    'CovariancePattern', 'Diagonal');

disp(glmm);

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
% 2. PLOTTING
%    All plotting loops iterate over the 'common' lists and look up
%    data by name, ensuring correct alignment.
% =========================================================================

ff = figure('color','w','Name','GLMM results');
col = get(0,'defaultAxesColorOrder');
ff.Units = 'centimeters';
ff.Position = [2 2.5 16 8];

ha = tight_subplot(1, 2, [.1 .1], [.15 .05], [.12 .03]);
for k=1:length(ha)
    ha(k).XTickLabelMode = 'auto';
    ha(k).YTickLabelMode = 'auto';
end

% --- PLOT 1: Betas and Sigmas ---
axes(ha(1)); hold on;
plot([0, n_betas_toplot+1], [0,0], 'k--');
plot_xticks = 1:n_betas_toplot;
plot_labels = cellfun(@(x) strrep(x, '_', '\_'), common_betanames, 'UniformOutput', false);
plot_labels = cellfun(@(x) strrep(x, '\_\_', '*'), plot_labels, 'UniformOutput', false);

for i = 1:n_betas_toplot
    pred_name = common_betanames{i};
    
    % Find the index for this predictor in this specific model's list
    beta_idx = find(strcmp(model_results.betanames, pred_name));
    beta_val = model_results.betas(beta_idx);

    % Get the Standard Error from the Coefficients table
    se_val = model_results.SE(beta_idx);

    % Get the sigma from the map. Default to 0 if no random slope exists.
    if isKey(model_results.sigma_map, pred_name)
        sigma_val = model_results.sigma_map(pred_name);
    else
        sigma_val = 0;
    end

    % Combine the two sources of variance
    total_sigma = sqrt(se_val.^2 + sigma_val.^2);

    offset = 0;

    % Plot point and error bar
    plot(i + offset, beta_val, '.', 'MarkerSize', 17, 'Color', col(1,:));
    if total_sigma > 0
        errbar(i + offset, beta_val, total_sigma, 'color', col(1,:));
        errbar(i + offset +.05, beta_val, sigma_val, 'color', col(1,:),'linestyle',':');
    end
end
ylabel('Betas'); xlim([0, n_betas_toplot+1]); view([90 90]);
ha(1).XTick = plot_xticks; ha(1).XTickLabel = plot_labels; ha(1).XGrid = 'on';

% --- PLOT 2: T-Statistics ---
axes(ha(2)); hold on;
plot([0, n_betas_toplot+1], [0,0], 'k--');
plot([0, n_betas_toplot+1], [2,2], 'k:'); plot([0, n_betas_toplot+1], [-2,-2], 'k:');

for i = 1:n_betas_toplot
    pred_name = common_betanames{i};
    tStat_idx = find(strcmp(model_results.betanames, pred_name));
    tStat_val = model_results.tStats(tStat_idx);
    % Plotting offset
    offset = 0;
    plot(i + offset, tStat_val, '.', 'MarkerSize', 17, 'Color', col(1,:));
end
ylabel('t-values'); xlim([0, n_betas_toplot+1]); view([90 90]);
ha(2).XTick = plot_xticks; ha(2).XTickLabel = plot_labels; ha(2).XGrid = 'on';
