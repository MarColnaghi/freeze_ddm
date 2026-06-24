clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion/proximity', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 150);
bouts_proc = bouts_proc(bouts_proc.durations < 630, :);

% Extract the pre-freezing social motion
inizio = 0;
fine = 630;
delay_after = 60;

distance_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine], 'cache', 'mindist_cache');
distance_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'align', 'offset', 'delay_after', delay_after, 'cache', 'mindist_cache');
distance_freeze_inv = 1./distance_freeze;

sm_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine]);
sm_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'align', 'offset', 'delay_after', delay_after);

ss_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine], 'cache', 'sumspeed_cache');
ss_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'align', 'offset', 'delay_after', delay_after, 'cache', 'sumspeed_cache');

sa_ = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine], 'cache', 'sumangle_cache');
sa_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'align', 'offset', 'delay_after', delay_after, 'cache', 'sumangle_cache');


%% Setup Datasets and Metadata
datasets = {distance_freeze, ss_freeze, sm_freeze, sa_freeze, distance_freeze_inv};
datasets_str = {'mindist', 'ss', 'sm', 'sa', 'mindist-1'};

y_labels = {'Min Distance', 'Social Speed', 'Social Motion', 'Angular Size', 'Min Distance ^-1'};
y_maxs = {700, 90, 25, 2.5, 0.035};

hex_colors = {'#6E44FF', '#B892FF', '#FFC2E2', '#EF7A85'};

% Time and Windowing Parameters
window_size = 3;
windows = [120, 60, 30, 0];
xlims = [-180 60];

for d = 1:numel(datasets)
    curr_data = datasets{d};
    
    % Initialize Figure
    fh = figure('color','w','Position',[100 100 1100 300]); 
    tlo = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % --- TILE 1 & 2: Mean Plot with Horizontal Distributions ---
    nexttile(1, [1 2])
    hold on
    t = (-size(curr_data, 2) + 1 : 1 : 0) + delay_after;
    
    % 1. Calculate Mean and SEM
    avg = mean(curr_data, 1, 'omitnan');
    sem = std(curr_data, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(curr_data), 1));
    
    % 2. Dynamic Y-Axis Limits (for breathing room)
    y_min = min(avg) - 0.1;
    y_max = max(avg) + 0.1;
    padding = (y_max - y_min) * 0;
    if padding == 0, padding = 1; end

    y_max = y_maxs{d};
    curr_ylim = [0, y_max];
    
    % 3. Define Histogram Edges for this specific dataset
    edges = linspace(curr_ylim(1), curr_ylim(2), 75);
    centers = edges(1:end-1) + diff(edges)/2;
    
    total_size_y = diff(curr_ylim);
    % 4. Plot Horizontal Window Distributions
    for i = 1:length(windows)
        w = windows(i);
        x0 = -w;
        
        % Extract window data
        idx = size(curr_data, 2) - delay_after - w;
        hvals = mean(curr_data(:, idx - window_size : idx + window_size), 2, 'omitnan');
        
        counts = histcounts(hvals, edges);
        counts = (counts / max(counts)) * 18; % Scale width relative to x-axis
        
        % Plot Vertical "Violin" Patches
        patch([x0 + counts, fliplr(x0*ones(size(counts)))], ...
              [centers, fliplr(centers)], ...
              hex2rgb(hex_colors{i}), 'FaceAlpha', 0.4, 'EdgeColor', 'k', 'LineWidth', 0.5);
          
        % Plot reference "marker" on the x-axis floor
        fill([x0-window_size x0-window_size x0+window_size x0+window_size], ...
             [curr_ylim(1) curr_ylim(1)-0.035*total_size_y curr_ylim(1)-0.035*total_size_y curr_ylim(1)], ...
             hex2rgb(hex_colors{i}), 'Clipping', 'off', 'EdgeColor', 'none');
    end
    
    % 5. Plot Main Line and SEM Shading
    fill([t, fliplr(t)], [avg+sem, fliplr(avg-sem)], 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    plot(t, avg, 'LineWidth', 2.5, 'Color', 'k')
    
    % Formatting Tile 1/2
    xline(0, 'k--', 'LineWidth', 1.5);
    ylabel(y_labels{d});
    xlabel('Frames relative to Freeze');
    apply_generic(gca, 'font_size', 18, 'tick_length', 0.015, 'line_width', 2, ...
        'xticks', sort([-180 -windows 60]), 'ylim', curr_ylim, 'xlim', xlims);
    set(gca, 'Layer', 'top')

    % --- TILE 3: PDF Comparison ---
    nexttile
    hold on
    for i = 1:length(windows)
        w = windows(i);
        idx = size(curr_data, 2) - delay_after - w;
        hvals = mean(curr_data(:, idx - window_size : idx + window_size), 2, 'omitnan');
        
        histogram(hvals, edges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', ...
            'EdgeColor', hex_colors{i}, 'LineWidth', 2);
    end
    
    view([90 -90]) % Flip this to match the orientation of the horizontal patches
    ylabel('Density');
    apply_generic(gca, 'font_size', 18, 'tick_length', 0.015, 'line_width', 2, 'xlim', curr_ylim);
    
    % Export
    save_name = sprintf('analysis_%s', datasets_str{d});
    exporter(fh, paths, [save_name, '.pdf']);
    exporter(fh, paths, [save_name, '.png']);

    fprintf('Exported %s\n', save_name);
end
