%% Bout-level Predictor Correlation Matrix (Pooled Across Flies)
%
% This script computes a correlation matrix across immobility bouts,
% treating each bout as an independent observation.
%
% Fly identity is intentionally ignored.
%
% Purpose:
%   - Diagnose collinearity between candidate predictors
%   - Identify redundant or bookkeeping-derived variables
%   - Sanity-check feature construction before model inclusion
%
% This analysis is NOT:
%   - A fly-level analysis
%   - A hierarchical or random-effects correlation
%   - Intended for inference about individual differences
%

disp('Generating bout-level predictor correlation heatmap...');

% =========================================================================
% 1. CONFIGURATION
% =========================================================================

id_code = 'imm3_mob3_pc4';
paths_out = path_generator('folder', 'experimental_vars/freeze_durations', 'bouts_id', id_code, 'imfirst', false);

thresholds = define_thresholds('min_dur', 30);
contact_threshold = 70;

% =========================================================================
% 2. LOAD AND PREPROCESS BOUT DATA
% =========================================================================

bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);

% The same freezes as the freeze duration analyses (fd_*): loom-evoked, first
% loom and already-frozen starts left out, contact-censored
bouts_proc = loom_evoked_freezes(bouts, 'already_frozen', false, 'min_dur', thresholds.min_dur, 'contact_threshold', contact_threshold);

% late_onset is 1 for the freezes whose onset is after the split latency of
% their loom speed, the dashed line in le_threshold_plot, in frames after the loom
split_latency = [mean(thresholds.le_window_sl) - 1.5, mean(thresholds.le_window_fl) - 2.5];
split_latency_row = reshape(split_latency(1 + (bouts_proc.sloom == 50)), [], 1); % column, indexing a row with a column gives a row
bouts_proc.late_onset = double(bouts_proc.onsets_loomaligned > split_latency_row);

T = bouts_proc;

% =========================================================================
% 3. EXTRACT BOUT-LEVEL PREDICTOR MATRIX
% =========================================================================
%
% Rows   = individual immobility bouts
% Columns = predictor variables
%
% Data are pooled across flies and loom presentations.

predictor_names = { 'avg_fs_loom',...
    'avg_fs',...
    'avg_fs_1s', ...
    'avg_pc',...
    'ls', ...
    'moving_flies', ...
    'avg_sm', ...
    'avg_ss', ...
    'time_since_last', ...
    'ln', ...
    'n_generated_freezes', ...
    'cum_freeze_time', ...
    'avg_history_dur', ...
    'late_onset'};

matrix = T{:, predictor_names};

% Bouts with a NaN in any predictor are dropped, so every pair of predictors is
% correlated over the same bouts
matrix = matrix(~any(isnan(matrix), 2), :);
fprintf('%d of %d bouts kept\n', size(matrix, 1), height(T));


% =========================================================================
% 4. COMPUTE BOUT-LEVEL CORRELATION MATRIX
% =========================================================================
%
% Correlations are computed across bouts, not across flies.
% This is intended for multicollinearity diagnostics only.

corr_matrix = corr(matrix, 'Rows', 'complete', 'Type', 'Pearson');

% Compute distance based on correlation
dist_metric = 1 - abs(corr_matrix);
link = linkage(squareform(dist_metric), 'ward');
leaf_order = optimalleaforder(link, squareform(dist_metric));

% Reorder everything
corr_matrix = corr_matrix(leaf_order, leaf_order);
predictor_names = predictor_names(leaf_order);

% =========================================================================
% 5. VISUALIZE CORRELATION MATRIX
% =========================================================================

fh = figure('Position', [100 100 800 800], 'Color', 'w');
t = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'loose');

imagesc(corr_matrix, [-1 1])
axis square
hold on
yticks(1:length(predictor_names))
xticks(1:length(predictor_names))

xticklabels(strrep(predictor_names,'_','\_'))
yticklabels(strrep(predictor_names,'_','\_'))
apply_generic(gca, 'box', 'on', 'font_size', 16, 'line_width', 1.5)

cmap = cbrewer2('div', 'RdBu', 60, 'pchip');
colormap(flipud(cmap))   % blue = negative, red = positive (usual convention)
cb = colorbar;
cb.Ticks = -1:0.5:1;
cb.TickDirection = 'out';
cb.LineWidth = 1.5;
cb.FontSize = 16;
cb.Label.String = 'Pearson r';
cb.Label.FontSize = 18;

n = size(corr_matrix, 1);

for i = 1:n
    for j = 1:n
        val = corr_matrix(i,j);

        % Skip NaNs just in case
        if isnan(val), continue, end

        % Choose text color for contrast
        if abs(val) > 0.6
            txt_col = 'w';
        else
            txt_col = 'k';
        end

        text(j, i, sprintf('%.3f', val), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 12, ...
            'Color', txt_col);
    end
end


exporter(fh, paths_out, 'corr_matrix.pdf')
