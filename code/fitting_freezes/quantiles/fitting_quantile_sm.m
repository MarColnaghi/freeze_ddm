% clearvars

% Importing the data
col = cmapper();
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_freezes/le/quantiles', 'bouts_id', id_code);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
bouts = bouts_formatting(bouts, thresholds);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 1:20, 'min_dur', 30);

% Extract the pre-freezing social motion
sm_pre = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [-25 5]);
bouts_proc.smp = mean(sm_pre, 2);

% Censoring and Truncation point for the fitting
points.censoring = 10.5;
points.truncation = 0;

% Distributions of the control variables
fh = figure('Position', [100 100 500 400], 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'loose')

nexttile
histogram(bouts_proc.fs, 0:0.02:3, 'FaceColor', col.vars.fs(round(end/2), :), 'EdgeColor', 'none')
quant.fs = prctile(bouts_proc.fs, [0, 25, 50, 75, 100]); quant.fs(1) = 0; quant.fs(end) = 2;
% quant.fs = [0 0.45 0.7 1.0 2.2];
xline(quant.fs);
apply_generic(gca, 'xlim', [0 2])

nexttile
histogram(bouts_proc.sm, 0:0.02:3, 'FaceColor', col.vars.sm(round(end/2), :), 'EdgeColor', 'none');
quant.sm = prctile(bouts_proc.sm, [0, 25, 50, 75, 100]); quant.sm(1) = 0; quant.sm(end) = 1.2;
% quant.sm = [0, .25, 0.45, 0.8, 1.85];
xline(quant.sm);
apply_generic(gca, 'xlim', [0 2])

quant.nloom_norm = [-0.05 0.55 1.05 1.55 2.05];

% Let's move onto the fitting procedure

model_2_fit = 'dddm0';
model_func = str2func(strcat('model_', model_2_fit));
model_struct = model_func();

n_quantiles = 4;
p = linspace(0, 100, n_quantiles + 1);

quant.sm = prctile(bouts_proc.sm, p);
n_looms = length(unique(bouts_proc.ls));
n_params = length(fieldnames(model_struct));

estimates = nan(n_quantiles, n_looms, n_params);
x = nan(n_quantiles, n_looms);

vars = {'sm'};

paths = path_generator('folder', fullfile('fitting_freezes','le','quantiles', model_2_fit), 'bouts_id', id_code);
%create_output_dirs(paths)

% --- Initialization ---
parameters = struct();
parameters.means = nan(n_quantiles, numel(unique(bouts_proc.ls)), length(fieldnames(model_struct)));
parameters.std = nan(n_quantiles, numel(unique(bouts_proc.ls)), length(fieldnames(model_struct)));
parameters.medians = nan(n_quantiles, numel(unique(bouts_proc.ls)), length(fieldnames(model_struct)));
parameters.IQR = nan(n_quantiles, numel(unique(bouts_proc.ls)), length(fieldnames(model_struct)), 2);

fh_main = figure('color','w','Position',[100, 100, 700, 320]);
tlo = tiledlayout(1, size(parameters.means, 3), 'TileSpacing', 'loose', 'Padding', 'compact');
x_vals = 1:n_quantiles;

for p = 1:size(parameters.means, 3)
    nexttile
    hold on; axis square;
    % Initialize with NaNs so the plots exist but are empty
    hSlow(p) = scatter(x_vals - 0.1, parameters.means(:, 1, p), 50, 'r', 'filled');
    hFast(p) = scatter(x_vals + 0.1, parameters.means(:, 2, p), 50, 'b', 'filled');
    
    hErrSlow(p) = errorbar(x_vals - 0.1, parameters.means(:,1,p), parameters.std(:,1,p), ...
        'o', 'Color', 'r', 'LineWidth', 1.2);
        
    hErrFast(p) = errorbar(x_vals + 0.1, parameters.means(:,2,p), parameters.std(:,2,p), ...
        'o', 'Color', 'b', 'LineWidth', 1.2);

    xlim([0.5, n_quantiles + 0.5]);
    xticks(x_vals);
    title(sprintf('Param %d', p)); % Better for tracking
    apply_generic(gca)
end

% --- Processing Loop ---
for idx_vars = vars
    vr = idx_vars{1};
    for idx_quantiles = 1:n_quantiles
        % Use a logic-based approach for ls (assuming 0 and 1)
        ls_values = unique(bouts_proc.ls)'; % e.g., [0, 1]
        
        for idx_ls = ls_values
            % 1. Data Processing & Fitting
            bouts_quant = quantilizer_v2(bouts_proc, ...
                'indexed_quantile', struct('sm', idx_quantiles, 'ls', idx_ls, 'fs', 1), ...
                'total_quantiles', struct('sm', n_quantiles, 'fs', 1));
            
            [model_results, ~] = run_fitting_newer(bouts_quant, points, model_2_fit, paths, ...
                'export', false, 'bads_display', 'off', 'n_bads', 3);
            
            % 2. Update the Parameter Matrix
            % Assuming idx_ls is 0-based, we use idx_ls + 1 for indexing
            col = idx_ls + 1; 
            parameters.means(idx_quantiles, col, :) = model_results.array_mean;
            parameters.std(idx_quantiles, col, :) = model_results.array_std;
 
            parameters.medians(idx_quantiles, col, :) = model_results.array_median;
            parameters.IQR(idx_quantiles, col, :, :) = model_results.array_IQR';

            % 3. Dynamic Plot Update
            for p = 1:size(parameters.means, 3)
                % Update Slow Loom (Column 1)
                set(hErrSlow(p), 'YData', parameters.means(:, 1, p));
                set(hErrSlow(p), 'YNegativeDelta', parameters.std(:, 1, p));
                set(hErrSlow(p), 'YPositiveDelta', parameters.std(:, 1, p));

                set(hErrSlow(p), 'YData', parameters.medians(:, 1, p));
                set(hErrSlow(p), 'YNegativeDelta', parameters.medians(:, 1, p) - parameters.IQR(:, 1, p, 1));
                set(hErrSlow(p), 'YPositiveDelta', parameters.IQR(:, 1, p, 2) - parameters.medians(:, 1, p));

                % Update Fast Loom (Column 2)
                set(hErrFast(p), 'YData', parameters.means(:, 2, p));
                set(hErrFast(p), 'YNegativeDelta', parameters.std(:, 2, p));
                set(hErrFast(p), 'YPositiveDelta', parameters.std(:, 2, p));

                set(hErrFast(p), 'YData', parameters.medians(:, 2, p));
                set(hErrFast(p), 'YNegativeDelta', parameters.medians(:, 2, p) - parameters.IQR(:, 2, p, 1));
                set(hErrFast(p), 'YPositiveDelta', parameters.IQR(:, 2, p, 2) - parameters.medians(:, 2, p));

            end
            drawnow % Force MATLAB to render the change
        end
    end
end