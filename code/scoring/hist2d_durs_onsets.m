%% Script Name:
% hist2d_durs_onsets
% last update: 28.07
%
% hist2d_durs_onsets('bouts', bouts, 'thresholds', thresholds, ...)
%
% Name/value parameters (all optional):
%   bouts            bout table. If empty, the dataset in paths.dataset is loaded
%   thresholds       thresholds struct                       (default: define_thresholds)
%   time_window      [pre post] in frames, pre positive      (default: thresholds.time_window)
%   type             'loom' | 'bsl'                          (default: 'loom')
%   str              'immobility' | 'mobility'               (default: 'immobility')
%   plot_style       'tile' | 'bar3'                         (default: 'tile')
%   split_by         bouts column to put one tile per level  (default: 'moving_flies')
%   zoom             logical                                 (default: false)
%   paths            paths struct, only needed for exporting (default: path_generator)
%   export           logical                                 (default: false)
%   genotype         genotype filter, auto-load only         (default: 1)
%   clamp            censor durations at the window end      (default: false)
%   bouts_with_loom  logical filter, or 'all' for no filter  (default: 'all')
%   censoring        bouts column flagging right-censored bouts, typically
%                    'is_censored' from impose_contact_threshold, which covers
%                    both collisions and the observation window. Adds a second
%                    row of tiles                            (default: '', off)
%   cens_style       'split'    top row uncensored, bottom row censored
%                    'fraction' top row everything, bottom row the censored
%                               fraction of each bin         (default: 'split')
%   dur_var          bouts column on the duration axis. Defaults to
%                    'ending_time' when censoring is on, since a censored bout
%                    is only observed up to that point, and 'ending_time'
%                    equals 'durations' for every uncensored bout
%                                                       (default: 'durations')
%   share_scale      put every tile, both rows included, on one colour and z
%                    scale                                   (default: true)
%   min_count        bins holding fewer bouts than this are left blank, in
%                    cens_style 'fraction' only              (default: 5)
%   cens_bin_size    bin width in frames for the fraction tiles. The count
%                    grid is far too fine to estimate a ratio on, so these
%                    bins are coarser                        (default: 15)

function hist2d_durs_onsets(varargin)

opt = inputParser;
addParameter(opt, 'bouts', []);
addParameter(opt, 'thresholds', []);
addParameter(opt, 'time_window', []);
addParameter(opt, 'type', 'loom', @(x) any(strcmp(x, {'loom', 'bsl'})));
addParameter(opt, 'str', 'immobility', @(x) any(strcmp(x, {'immobility', 'mobility'})));
addParameter(opt, 'plot_style', 'tile', @(x) any(strcmp(x, {'tile', 'bar3'})));
addParameter(opt, 'split_by', 'moving_flies');
addParameter(opt, 'zoom', false);
addParameter(opt, 'paths', []);
addParameter(opt, 'export', false);
addParameter(opt, 'genotype', 1);
addParameter(opt, 'clamp', false);
addParameter(opt, 'bouts_with_loom', 'all');
addParameter(opt, 'censoring', '');
addParameter(opt, 'cens_style', 'split', @(x) any(strcmp(x, {'split', 'fraction'})));
addParameter(opt, 'dur_var', 'durations');
addParameter(opt, 'share_scale', true);
addParameter(opt, 'min_count', 5);
addParameter(opt, 'cens_bin_size', 15);

parse(opt, varargin{:});

bouts = opt.Results.bouts;
thresholds = opt.Results.thresholds;
type = opt.Results.type;
str = opt.Results.str;
plot_style = opt.Results.plot_style;
split_by = opt.Results.split_by;
zoom_flag = opt.Results.zoom;
paths = opt.Results.paths;
export = opt.Results.export;
clamp = opt.Results.clamp;
bouts_with_loom = opt.Results.bouts_with_loom;
censoring = opt.Results.censoring;
cens_style = opt.Results.cens_style;
share_scale = opt.Results.share_scale;
min_count = opt.Results.min_count;
cens_bin_size = opt.Results.cens_bin_size;

% Load Colors
extra.quantiles = 5;
col = cmapper([], extra.quantiles);

% Define thresholds
if isempty(thresholds)
    thresholds = define_thresholds;
end

time_window = opt.Results.time_window;
if isempty(time_window)
    time_window = thresholds.time_window;
end

% Load Paths & Files
if isempty(paths) && (export || isempty(bouts))
    paths = path_generator('folder', 'validation_scoring/durs_onsets');
end

if isempty(bouts)
    disp('no bouts provided, will load the dataset in the dataset folder')

    loaded = load(fullfile(paths.dataset, 'bouts.mat'));
    bouts = loaded.bouts;
    bouts = bouts(bouts.genotype == opt.Results.genotype, :);
    bouts = bouts_formatting(bouts, thresholds);

else
    disp('bouts provided')
end

% The censoring flag is not part of the bouts table as it comes out of
% load_flies_new, so fail here rather than half way through the layout
if ~isempty(censoring) && ~ismember(censoring, bouts.Properties.VariableNames)
    error('hist2d_durs_onsets:missingCensorColumn', ...
        ['Column ''%s'' is not in the bouts table. Add it with ' ...
        'impose_contact_threshold, e.g. bouts = impose_contact_threshold(bouts);'], censoring);
end

% clamp truncates every duration at the window end, which is a censoring
% scheme of its own. Running it underneath another one hides which boundary
% actually stopped each bout
if clamp && ~isempty(censoring)
    error('hist2d_durs_onsets:clampWithCensoring', ...
        'clamp and censoring (''%s'') are two censoring schemes. Pick one.', censoring);
end

% A censored bout is only observed up to ending_time, so that is where it
% belongs on the duration axis. ending_time already equals durations for every
% bout that was never censored, so one column serves both rows
dur_var = opt.Results.dur_var;
if ~isempty(censoring) && ismember('dur_var', opt.UsingDefaults) && ...
        ismember('ending_time', bouts.Properties.VariableNames)
    dur_var = 'ending_time';
end

if ~ismember(dur_var, bouts.Properties.VariableNames)
    error('hist2d_durs_onsets:missingDurationColumn', ...
        'Column ''%s'' is not in the bouts table.', dur_var);
end

switch plot_style
    case 'bar3'
        bin_size = 1;
    case 'tile'
        bin_size = 2;
end

% Filters that do not depend on the split variable, applied once. They are
% all row-wise, so pulling them out of the loop leaves each tile identical
% to what the per-condition filtering produced
if clamp == true
    bouts.durations(bouts.durations + bouts.onsets_loomwin >= time_window(2)) = time_window(2) - bouts.onsets_loomwin(bouts.durations + bouts.onsets_loomwin >= time_window(2));
end

if islogical(bouts_with_loom)
    bouts = bouts(bouts.bout_with_loom == bouts_with_loom, :);
end

if strcmp(type, 'loom')
    bouts = bouts(bouts.period == 1, :);

elseif strcmp(type, 'bsl')
    bouts = bouts(bouts.period == 0, :);
end

if strcmp(str, 'immobility')
    bouts = bouts(bouts.type == 1, :);

elseif strcmp(str, 'mobility')
    bouts = bouts(bouts.type == 0, :);
end

conditions = unique(bouts.(split_by))';
n_cond = numel(conditions);

durs_all = bouts.(dur_var);

% One row of tiles per censoring status. Bouts whose status is unknown are
% dropped from both rows rather than quietly counted as uncensored
if isempty(censoring)
    row_mask = {true(height(bouts), 1)};
    row_name = {''};

elseif strcmp(cens_style, 'split')
    cens_col = double(bouts.(censoring));
    known = ~isnan(cens_col);
    row_mask = {known & cens_col ~= 1, known & cens_col == 1};
    row_name = {'uncensored', 'censored'};

else
    % The fraction row is not a histogram, so the single histogram row
    % carries every bout, censored or not
    row_mask = {true(height(bouts), 1)};
    row_name = {'all bouts'};
end

n_hist_rows = numel(row_mask);           % rows drawn as histogram2
n_rows = 1 + double(~isempty(censoring)); % rows of tiles in the layout

% Axis limits are needed up front, so that only the bins that land
% inside them get handed to the renderer
if zoom_flag
    xl = [-1 81];
    yl = [-10 61];
else
    xl = [-10 900];
    yl = [-time_window(1) - 10, time_window(2) + 25];
end

% One bin grid for every tile, so that the tiles are directly comparable
x_edges_full = -0.5:bin_size:max(durs_all) + 0.5;
y_edges_full = -time_window(1) - 0.5:bin_size:time_window(2) + 0.5;

% Counts over the full grid of every tile, so that neither the cropping
% below nor a sparse tile can move the colour and z scaling
max_count = zeros(1, n_hist_rows);
for r = 1:n_hist_rows
    for idx_cond = conditions
        mask = row_mask{r} & bouts.(split_by) == idx_cond;
        max_count(r) = max(max_count(r), max(histcounts2(durs_all(mask), bouts.onsets_loomwin(mask), x_edges_full, y_edges_full), [], 'all'));
    end
end

% Censored bouts are a minority, so a shared scale makes the bottom row read
% as faint. That is the honest picture of how many there are, but it hides
% the shape, hence the switch
if share_scale
    max_count(:) = max(max_count);
end

max_count = max(max_count, 1); % keeps clim valid when a selection is empty

% Keep only the bins inside the axis limits. The bin boundaries are
% untouched, so every visible bin holds exactly the same count as before
x_edges = x_edges_full(x_edges_full >= xl(1) - bin_size & x_edges_full <= xl(2) + bin_size);
y_edges = y_edges_full(y_edges_full >= yl(1) - bin_size & y_edges_full <= yl(2) + bin_size);

% The fraction tiles need their own, coarser grid: a count of one is a
% perfectly good count, but a ratio out of one bout is only ever 0 or 1.
% Every fraction bin is a whole number of count bins, so the rows stay
% aligned, they just no longer match one for one. A trailing group that does
% not fill k_cens count bins is dropped, which keeps the bins uniform
k_cens = max(1, round(cens_bin_size / bin_size));
x_edges_cens = x_edges(1:k_cens:end);
y_edges_cens = y_edges(1:k_cens:end);

% Create Figure
fh = figure('color','w','Position',[100, 100, 350 * n_cond, 400 * n_rows + 350]);
tl = tiledlayout(fh, n_rows, n_cond, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = gobjects(n_rows, n_cond);

for r = 1:n_hist_rows
    i = 0;
    for idx_cond = conditions
        i = i + 1;
        mask = row_mask{r} & bouts.(split_by) == idx_cond;
        % mask = mask & bouts.nloom == 16;

        ax(r, i) = nexttile((r - 1) * n_cond + i);
        hold on

        h2h = histogram2(durs_all(mask), bouts.onsets_loomwin(mask), x_edges, y_edges, 'DisplayStyle', plot_style, 'FaceColor', 'flat', 'EdgeColor', 'none');
        % If you want marginals
        % h2h = histogram2(durs, 604 * ones(length(durs), 1), -0.5:1:601.5, 603:605, 'DisplayStyle', 'bar3', 'FaceColor', 'k', 'Normalization', 'pdf', 'EdgeColor','none');
        % h2h = histogram2(-1 * ones(length(onsets), 1), onsets, -2:0,  -0.5:1:601.5, 'DisplayStyle', 'bar3', 'FaceColor', 'k', 'Normalization', 'pdf', 'EdgeColor','none');

        % Set Axis Labels
        if strcmp(dur_var, 'durations')
            xlabel('Bout Duration (frames)');
        else
            % Says outright that a censored bar sits at the point observation
            % stopped, not at the length the bout would have reached
            xlabel(sprintf('Bout Duration (frames, %s)', strrep(dur_var, '_', ' ')));
        end
        ylabel('Bout Latency (frames)');
        zlabel('Count');

        if isempty(row_name{r})
            title(sprintf('%s = %g', strrep(split_by, '_', ' '), idx_cond));
        else
            title(sprintf('%s = %g, %s (n = %d)', strrep(split_by, '_', ' '), idx_cond, row_name{r}, sum(mask)));
        end

        % Set Axis Limits
        xlim(xl);
        ylim(yl);
        zlim([-0.2, max_count(r) + 10]);

        % Set Colormap
        colorcet('L08');
        clim([0 max_count(r)])

        % Apply Generic Changes
        apply_generic(ax(r, i), 'font_size', 20)
        set(ax(r, i),'TickLength',[0.02, 0.02])
        grid on

        % Keep drag to rotate and scroll to zoom, but drop the data tip, which
        % hit-tests against every face of the histogram on mouse move.
        % Pan is left out on purpose: it binds to drag as well and would
        % conflict with the rotation
        ax(r, i).Interactions = [rotateInteraction zoomInteraction];

        switch plot_style
            case 'bar3'
                thresholds.sp_window = [];
              %  ax.XLabel.Rotation = 300;
              %  ax.XLabel.Position(1:2) = [250 -225];
               % ax.YLabel.Rotation = 7;
              %  ax.YLabel.Position(1:2) = [675 210];
                if ~zoom_flag
                    view(80, 50)
                else
                    view(80, 30)

                end
                % Set Colormap
                colorcet('L08');
                clim([0 max_count(r)])

            case 'tile'
                bin_size = 2;
                view(90, 90)
                % Set Colormap
                colorcet('L02');
                clim([0 max_count(r)])

        end

        % Set Axis Ticks
        if ~zoom_flag
            xticks(0:100:900)
        else
            xticks([0 60])
        end
        ticks = sort([-time_window(1) 0 30 100 200 300 400 500 time_window(2)]);
        yticks(ticks)
        ytlh = yticklabels; ytlh{2} = 'onset'; ytlh{3} = 'offset'; yticklabels(ytlh); ytickangle(360 - 60);

        % Create Plane of Loom Presentation
        x = meshgrid(-10:905:900);
        y = meshgrid(0:30:30)';
        z = [-0.1 -0.1; -0.1 -0.1];

        zmin = ax(r, i).ZLim(1);
        zmax = ax(r, i).ZLim(2);
        zrange = zmax - zmin;

        z_label_pos = zmax + 0.05 * zrange;  % 5% above the top of the plot

        lp = gobjects(0); % reset per tile, so nothing leaks into the next one

        switch type
            case 'loom'
                % Loom speeds are pooled, so the plane is drawn in a neutral
                % colour unless the tiles happen to be split by speed
                if strcmp(split_by, 'sloom') && idx_cond == 25
                    plane_col = col.vars.sloom(3,:); plane_txt = {'Slow', 'Loom'};
                elseif strcmp(split_by, 'sloom') && idx_cond == 50
                    plane_col = col.vars.sloom(5,:); plane_txt = {'Fast', 'Loom'};
                else
                    plane_col = col.loom; plane_txt = {'Loom'};
                end

                lp = surf(x, y, z, 'FaceColor', plane_col, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'LineWidth', 1.5);
                surf(z, y, x, 'FaceColor', plane_col, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'LineWidth', 1.5);
                text(0, 15, z_label_pos, plane_txt, 'FontSize', 18, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Rotation', ax(r, i).YLabel.Rotation)

            case 'bsl'
                text(0, 15, z_label_pos, {'Surrogate', 'Loom'}, 'FontSize', 18, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Rotation', ax(r, i).YLabel.Rotation)

        end

        if ~isempty(thresholds.sp_window)
            set(lp, 'Visible', 'off');
            for idx_window = 1:size(thresholds.sp_window, 1)
                if thresholds.sp_window(idx_window, 1) == 45
                    col.areas(end - idx_window, :) = [87 167 115]/255;
                end
                surf(meshgrid(-5:sum(time_window):sum(time_window)), meshgrid(thresholds.sp_window(idx_window, 1):30:thresholds.sp_window(idx_window, 2))',...
                    [-0.01 -0.01; -0.01 -0.01], 'FaceColor', col.areas(end - idx_window, :), 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'LineWidth', 1);
            end
        end
        axis square

    end
end

% Censored fraction row. Collision censoring depends on where the other flies
% were, so unlike the window and next loom boundaries it cannot be read off
% the two axes already plotted
if ~isempty(censoring) && strcmp(cens_style, 'fraction')

    cens_col = double(bouts.(censoring));

    i = 0;
    for idx_cond = conditions
        i = i + 1;

        % Bouts whose status is unknown leave both the numerator and the
        % denominator, so that they cannot drag the fraction towards zero
        cond_mask = bouts.(split_by) == idx_cond;
        known = cond_mask & ~isnan(cens_col);
        cens = known & cens_col == 1;

        n_known = histcounts2(durs_all(known), bouts.onsets_loomwin(known), x_edges_cens, y_edges_cens);
        n_cens = histcounts2(durs_all(cens), bouts.onsets_loomwin(cens), x_edges_cens, y_edges_cens);

        frac = n_cens ./ n_known;
        frac(n_known < min_count) = NaN; % a handful of bouts only ever reads as 0 or 1

        % Says whether a blank row means no collisions or bins too sparse
        fprintf('%s = %g: %d of %d fraction bins hold >= %d bouts, %.1f%% of bouts censored\n', ...
            split_by, idx_cond, sum(n_known(:) >= min_count), numel(n_known), min_count, ...
            100 * sum(cens) / max(sum(known), 1));

        ax(2, i) = nexttile(n_cond + i);

        % histogram2 maps colour to count, so the fraction is drawn as an
        % image instead. Blank bins stay transparent rather than reading as 0
        im = imagesc(x_edges_cens(1:end - 1) + k_cens * bin_size/2, y_edges_cens(1:end - 1) + k_cens * bin_size/2, frac');
        im.AlphaData = ~isnan(frac');

        % imagesc flips the y axis and pins the data aspect ratio, neither of
        % which suits an axis pair in frames. Blank bins fall through to the
        % white figure, since apply_generic sets the axes Color to none below
        set(ax(2, i), 'YDir', 'normal', 'DataAspectRatioMode', 'auto')
        hold on

        xlabel('Bout Duration (frames)');
        ylabel('Bout Latency (frames)');
        title(sprintf('censored fraction (n \\geq %d)', min_count));

        xlim(xl);
        ylim(yl);

        colorcet('L03');
        clim([0 1])

        % Loom onset and offset, the only reference the counts row shares
        yline(0, 'k-', 'LineWidth', 1);
        yline(30, 'k-', 'LineWidth', 1);

        apply_generic(ax(2, i), 'font_size', 20)
        set(ax(2, i), 'TickLength', [0.02, 0.02])

        ticks = sort([-time_window(1) 0 30 100 200 300 400 500 time_window(2)]);
        xticks(ax(2, i), get(ax(1, i), 'XTick'))
        yticks(ax(2, i), ticks)
        lbl = yticklabels(ax(2, i)); lbl{2} = 'onset'; lbl{3} = 'offset'; yticklabels(ax(2, i), lbl);
        ytickangle(ax(2, i), 360 - 60)

        % Matches the orientation of the counts row when plot_style is tile
        view(90, 90)
        axis square

        ax(2, i).Interactions = [zoomInteraction panInteraction];
    end
end

if share_scale && ~strcmp(cens_style, 'fraction')
    % One bar for the whole layout, since every tile is on the same scale
    cb = colorbar(ax(1, end));
    cb.Layout.Tile = 'east';
    cb.Location = 'eastoutside';
    cb.Label.String = 'Count';
    cb.LineWidth = 2;
    cb.FontSize = 20;

else
    % A bar per row, attached to the row rather than to the layout, since the
    % rows are not on the same scale
    for r = 1:n_rows
        cb = colorbar(ax(r, end));
        cb.LineWidth = 2;
        cb.FontSize = 20;

        if r > n_hist_rows
            cb.Label.String = sprintf('Fraction censored (%s)', strrep(censoring, '_', ' '));
        else
            cb.Label.String = 'Count';
        end
    end
end

% Export Figure

if export
    figure_title = sprintf('%s_%s_%s_by_%s_z%d', plot_style, str, type, split_by, zoom_flag);

    if ~isempty(censoring)
        figure_title = [figure_title, '_', censoring];
    end

    exporter(fh, paths, [figure_title, '.png'])
   % exporter(fh, paths, [figure_title, '.pdf'])
end