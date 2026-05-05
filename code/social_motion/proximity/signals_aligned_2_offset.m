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
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);
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

%%

% 1. Setup Data and Metadata
datasets = {sm_freeze, ss_freeze, sa_freeze,distance_freeze_inv };
datasets_str = {'sm', 'ss', 'sa', 'mindist'};
y_labels = {'Social Motion', 'Social Speed', 'Angular Size', 'Min Dist ^-1'};

col = cmapper;

% 2. Initialize Figure
fh = figure('color','w','Position',[100 100 700 900]); 
tlo = tiledlayout(numel(datasets), 1 , 'TileSpacing', 'compact', 'Padding', 'compact');

xlims = [-300 60];

% 3. The Loop
for i = 1:numel(datasets)
    nexttile
    hold on
    
    ax(i) = gca;
    current_data = datasets{i};
    
    % 1. Calculate Time Vector
    t = (-size(current_data, 2) + 1 : 1 : 0) + delay_after;
    
    % 2. Calculate Stats for SEM
    avg = mean(current_data, 1, 'omitnan');
    sem = std(current_data, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(current_data), 1));
    
    % 3. Plot Shaded SEM (using fill)
    % We concatenate t and the bounds to create a closed polygon
    fill_x = [t, fliplr(t)];
    fill_y = [avg + sem, fliplr(avg - sem)];
    
    fill(fill_x, fill_y, col.var.(datasets_str{i}), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.2);
    
    % 4. Plot the Mean Line
    plot(t, avg, 'LineWidth', 2, 'Color', col.var.(datasets_str{i}))
    

    % 6. Formatting
    apply_generic(gca, 'font_size', 20, 'tick_length', 0.015, 'line_width', 2, ...
        'xticks', xlims(1):60:xlims(2), 'xlim', xlims);

    current_ticks = ax(i).XAxis.TickValues;

    % Convert ticks to strings
    labels = arrayfun(@num2str, current_ticks, 'UniformOutput', false);

    % Replace the label where tick == 0
    labels(current_ticks == 0) = {'Offset'};

    % Apply the new labels
    xticklabels(ax, labels)

    ylabel(y_labels{i})

    xline(0, '--k');
end

linkaxes(ax, 'x')

exporter(fh, paths, 'signals_aligned_2_offset.pdf')
exporter(fh, paths, 'signals_aligned_2_offset.png')
% The Main Loop


