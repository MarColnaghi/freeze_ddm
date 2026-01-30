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
paths_out = path_generator('folder', 'experimental_vars', 'bouts_id', id_code, 'imfirst', false);

% =========================================================================
% 1. CONFIGURATION
% =========================================================================

id_code = 'imm3_mob3_pc4';

% Threshold definitions
thresholds = define_thresholds;
thresholds.le_window_fl = [5 40];
thresholds.le_window_sl = [15 50];

% Duration limits (seconds)
Tmax   = 10.5;
Ttrunc = 0.5;

% =========================================================================
% 2. LOAD AND PREPROCESS BOUT DATA
% =========================================================================

bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);

% Select immobility bouts during loom period, late window
bouts_le = data_parser_new( ...
    bouts, ...
    'period', 'loom', ...
    'window', 'le', ...
    'type', 'immobility', ...
    'nloom', 2:20, ...
    'min_dur', Ttrunc * 60);

T = bouts_le;

% =========================================================================
% 3. EXTRACT BOUT-LEVEL PREDICTOR MATRIX
% =========================================================================
%
% Rows   = individual immobility bouts
% Columns = predictor variables
%
% Data are pooled across flies and loom presentations.

predictor_names = { ...
    'avg_fs_1s', ...
    'sloom', ...
    'moving_flies', ...
    'avg_sm', ...
    'avg_ss', ...
    'time_since_last', ...
    'nloom', ...
    'n_generated_freezes', ...
    'cum_freeze_time', ...
    'avg_history_dur'};

matrix = T{:, predictor_names};

% Remove rows with any NaNs (pairwise deletion would obscure structure)
[nan_rows, ~] = find(isnan(matrix));
matrix(unique(nan_rows), :) = [];


% =========================================================================
% 4. COMPUTE BOUT-LEVEL CORRELATION MATRIX
% =========================================================================
%
% Correlations are computed across bouts, not across flies.
% This is intended for multicollinearity diagnostics only.

[corr_matrix, P] = corr(matrix, 'Rows', 'Pairwise', 'Type', 'Pearson');

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

% Pretty labels (escape underscores for heatmap)
labels = cellfun(@(x) strrep(x, '_', '\_'), ...
                 predictor_names, ...
                 'UniformOutput', false);

labels = cellfun(@(x) strrep(x, '\_\_', '*'), ...
                 labels, ...
                 'UniformOutput', false);

imagesc(corr_matrix, [-1 1])
axis square
hold on
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


disp('Bout-level correlation matrix plotted successfully.');
exporter(fh, paths_out, 'corr_matrix.pdf')