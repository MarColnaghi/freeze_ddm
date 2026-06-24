clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
% thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);
threshold = 80;
[bouts_proc, contact_mask, is_below_threshold] = impose_contact_threshold(bouts_proc, 'threshold', threshold);
censoring = bouts_proc.durations < 630;
bouts_no_contacts = bouts_proc(~contact_mask & censoring, :);
bouts_with_contacts = bouts_proc(contact_mask & censoring, :);

%%

inizio = 0;
fine = 630;

sm_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'norm_factor', 10, 'cache', 'motion_cache','window', [0 630]);
angle_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'norm_factor', 2, 'cache', 'sumangle_cache','window', [0 630]);

angle_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'norm_factor', 2, 'cache', 'sumangle_cache');
angle_freeze = angle_freeze(:, 1:631);

sm_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'norm_factor', 10, 'cache', 'motion_cache');
sm_freeze = sm_freeze(:, 1:631);

% sm_ili = sm_ili .* ~is_below_threshold;

bouts_proc = bouts_proc(~contact_mask, :);
sm_ili = sm_ili(~contact_mask, :);
angle_ili = angle_ili(~contact_mask, :);
sm_freeze = sm_freeze(~contact_mask, :);
angle_freeze = angle_freeze(~contact_mask, :);

n_quantiles = 4;
total_size = nan(n_quantiles, 1);

for idx_quantiles_sm = 1:n_quantiles
    
    % 1. Data extraction
    [~, mask] = quantilizer_v2(bouts_proc, 'total_quantiles', ...
        struct('sm', 4, 'fs', 1), 'indexed_quantile', ...
        struct('fs', 1, 'sm', idx_quantiles_sm));

    total_size(idx_quantiles_sm) = numel(mask);
    
end

boundaries = [0; cumsum(total_size)];

col = cmapper();
[~, sort_order] = sort(bouts_proc.sm);

% Visualize the Heatmap
fh = figure('color','w','Position',[100, 100, 600, 4000]);
tiledlayout(4, 1, 'Padding', 'compact', 'TileSpacing', 'tight')
nexttile(2, [3 1])
ax_distr(2) = gca;

hold on
x_axis = 1:size(sm_ili, 2);
y_axis = 1:size(sm_ili, 1);
h = imagesc(x_axis, y_axis, sm_ili(sort_order, :), [0 1.5]);
set(h, 'AlphaData', ~isnan(sm_ili(sort_order, :)));

scatter(bouts_proc.durations(sort_order), 1:height(bouts_proc), 2, '|', 'k')
apply_generic(gca, 'xlim', [0 630], 'no_y', true, 'ylim', [-100 size(sm_ili, 1) + 100], 'xtick', 0:120:600, 'font_size', 20, 'line_width', 2)
colormap(cbrewer2('Reds',[]));

%cb = colorbar(ax, 'Location', 'southoutside', 'FontSize', 18, 'LineWidth', 2);
%cb.Label.String = 'Social Motion';
xticklabels({'Onset', '2', '4', '6', '8', '10'});

% Set the interpreter to TeX so it recognizes the \newline command
ax_distr(2).TickLabelInterpreter = 'tex';
xtickangle(0);

xlabel('Time (s)');

% cb = colorbar(ax, 'Location', 'southoutside', 'FontSize', 18, 'LineWidth', 2);
% cb.Label.String = 'Social Motion';

for idx_cluster = 1:n_quantiles
    fill([ax_distr(2).XLim(1)-2,ax_distr(2).XLim(1)-2,ax_distr(2).XLim(1)-17,ax_distr(2).XLim(1)-17], [boundaries(idx_cluster) boundaries(idx_cluster + 1) boundaries(idx_cluster + 1)   boundaries(idx_cluster)], ...
        col.vars.sm(idx_cluster,:), 'EdgeColor','none', 'Clipping', 'off');
    text(mean([ax_distr(2).XLim(1)-2,ax_distr(2).XLim(1)-17]), mean([boundaries(idx_cluster); boundaries(idx_cluster + 1)]), num2str(idx_cluster),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    hold on

end

% Add horizontal lines to separate clusters
hold on;
for idx_cluster = 1:length(boundaries)
    yline(boundaries(idx_cluster), 'k--', 'LineWidth', 1.1);
end

nexttile(1)
colors = col.vars.sm; 
hold on
for idx_quantiles_sm = 1:n_quantiles
    
    % 1. Data extraction
    [~, mask] = quantilizer_v2(bouts_proc, 'total_quantiles', ...
        struct('sm', 4, 'fs', 1), 'indexed_quantile', ...
        struct('ls', 0, 'fs', 1, 'sm', idx_quantiles_sm));
    
    if ~any(mask); continue; end
    
    data = sm_ili(mask, :);
    
    % 2. Calculate Mean and SEM
    mu  = mean(data, 1, 'omitnan');
    sem = std(data, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data), 1));
    
    % 3. Prepare vectors for fill (removing NaNs for a clean patch)
    x = 1:length(mu);
    valid = ~isnan(mu) & ~isnan(sem); % Ensure no breaks in the polygon
    
    x_vec = x(valid);
    upper = mu(valid) + sem(valid);
    lower = mu(valid) - sem(valid);
    
    % Create the closed loop: [x, fliplr(x)] and [upper, fliplr(lower)]
    X_patch = [x_vec, fliplr(x_vec)];
    Y_patch = [upper, fliplr(lower)];
    
    % 4. Plot the Shaded Area
    fill(X_patch, Y_patch, colors(idx_quantiles_sm + 1, :), ...
        'FaceAlpha', 0.3, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off'); 
    
    % 5. Plot the Mean Line
    plot(x_vec, mu(valid), 'Color', colors(idx_quantiles_sm + 1, :), ...
        'LineWidth', 1);
    
end

% Formatting
ax_distr(1) = gca;
apply_generic(gca, 'xlim', [0 630], 'ylim', [0 1.5], 'xtick', 0:120:600, 'font_size', 20, 'line_width', 2, 'no_x', true)

ylabel({'Social' ,' Motion'});

linkaxes([ax_distr(:)], 'x')
xlim([0 630]);
xticklabels({'0', '2', '4', '6', '8', '10'});

exporter(fh, paths, 'clusters_by_sm_quantiles.pdf')
exporter(fh, paths, 'clusters_by_sm_quantiles.png')


increment = 60;
lags = 0:increment:600;
y_top = -10; 
y_bot = -100;
box_width = increment; % This matches your lag increment
ax = gca;

col = cmapper([], numel(lags));

axes(ax_distr(2))
clim([91 91.2])

for idx_window = 1:numel(lags)

    curr_lag = lags(idx_window);

    % Draw the rectangle for this specific window
    fill([curr_lag, curr_lag + increment, curr_lag + increment, curr_lag], [y_top, y_top, y_bot, y_bot], ...
        col.pca(idx_window, :), 'EdgeColor', 'none', 'Clipping', 'off');

    % Add the window index or lag value text
    text(median([curr_lag, curr_lag + increment]), mean([y_top, y_bot]), num2str(idx_window),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 10, 'Color', 'w', 'FontWeight', 'bold');
end

exporter(fh, paths, 'clusters_by_sm_quantiles_nocolor.pdf')
exporter(fh, paths, 'clusters_by_sm_quantiles_nocolor.png')

%%

fh = figure('color','w','Position',[100, 100, 600, 400]);
hold on


% 2. Add the window color bar at the bottom
% Adjust these values based on how much space you want below the axis
y_top = 0; 
y_bot = +0.035;
box_width = increment; % This matches your lag increment
ax = gca;

for idx_window = 1:numel(lags)

    curr_lag = lags(idx_window);
    duration_mask = bouts_proc.durations > curr_lag & bouts_proc.durations <= curr_lag + increment;
    sm_freeze_only = sm_freeze(duration_mask, :);

    plot(x, mean(sm_freeze_only, 'omitnan'), 'Color', col.pca(idx_window, :), 'LineWidth', 2)

    % Draw the rectangle for this specific window
    fill([curr_lag, curr_lag + increment, curr_lag + increment, curr_lag], [y_top, y_top, y_bot, y_bot], ...
        col.pca(idx_window, :), 'EdgeColor', 'none', 'Clipping', 'off');

    % Add the window index or lag value text
    text(median([curr_lag, curr_lag + increment]), mean([y_top, y_bot]), num2str(idx_window),...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 10, 'Color', 'w', 'FontWeight', 'bold');
end
apply_generic(ax, 'xlim', [0 630], 'xtick', 0:120:600, 'font_size', 20, 'line_width', 2, 'ylim', [0 0.8])
set(ax ,'Layer', 'Top')
ylabel('Social Motion')
%,

exporter(fh, paths, sprintf('clusters_by_breaktime_%d.pdf', threshold))
exporter(fh, paths, sprintf('clusters_by_breaktime_%d.png', threshold))

%%
% Process the data: set thresholds for sm, fs and ln. Set minimum duration.

fh = figure('color','w','Position', [100, 100, 700, 300]);
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for idx_threshold = [0, 70]

    bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);
    threshold = idx_threshold;
    [bouts_proc, contact_mask] = impose_contact_threshold(bouts_proc, 'threshold', threshold);
    censoring = bouts_proc.durations < 630;

    sm_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'norm_factor', 1, 'cache', 'sumangle_cache', 'align', 'offset');
    bouts_proc = bouts_proc(~contact_mask & censoring, :);
    sm_ili = sm_ili(~contact_mask & censoring, :);

    col = cmapper();

    nexttile;
    hold on

    colors = col.vars.sm;

    for idx_quantiles_sm = 1:4

        % 1. Data extraction
        [~, mask] = quantilizer_v2(bouts_proc, 'total_quantiles', ...
            struct('sm', 4, 'fs', 1), 'indexed_quantile', ...
            struct('fs', 1, 'sm', idx_quantiles_sm));

        if ~any(mask); continue; end

        data = sm_ili(mask, :);

        % 2. Calculate Mean and SEM
        mu  = mean(data, 1, 'omitnan');
        sem = std(data, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data), 1));

        % 3. Prepare vectors for fill (removing NaNs for a clean patch)
        x = -length(mu):1:-1;
        valid = ~isnan(mu) & ~isnan(sem); % Ensure no breaks in the polygon

        x_vec = x(valid);
        upper = mu(valid) + sem(valid);
        lower = mu(valid) - sem(valid);

        % Create the closed loop: [x, fliplr(x)] and [upper, fliplr(lower)]
        X_patch = [x_vec, fliplr(x_vec)];
        Y_patch = [upper, fliplr(lower)];

        % 4. Plot the Shaded Area
        fill(X_patch, Y_patch, colors(idx_quantiles_sm + 1, :), ...
            'FaceAlpha', 0.3, ...
            'EdgeColor', 'none', ...
            'HandleVisibility', 'off');

        % 5. Plot the Mean Line
        plot(x_vec, mu(valid), 'Color', colors(idx_quantiles_sm + 1, :), ...
            'LineWidth', 2.5);

    end

    % Formatting
    apply_generic(gca, 'font_size', 20, 'tick_length', 0.015, 'line_width', 2, 'ylim', [0 1.5], 'xtick', -600:120:0);
    xlim([-630 0]);
    xtickangle(0)
    ylabel('Summed Angle');
    xlabel('Time (frames)');
end

exporter(fh, paths, 'clusters_by_sm_angle_offset_aligned.pdf')
exporter(fh, paths, 'clusters_by_sm_angle_offset_aligned.png')