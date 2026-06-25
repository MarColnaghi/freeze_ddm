clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion/event_locked', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
%thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);

% Extract the pre-freezing social motion
inizio = 0;
fine = 630;
delay_after = 90;

distance_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine], 'cache', 'mindist_cache');
distance_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'align', 'offset', 'delay_after', delay_after, 'cache', 'mindist_cache');
distance_freeze_inv = 1./distance_freeze;

sm_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine]);
sm_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'align', 'offset', 'delay_after', delay_after);

ss_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine], 'cache', 'sumspeed_cache');
ss_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'align', 'offset', 'delay_after', delay_after, 'cache', 'sumspeed_cache');

sa_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine], 'cache', 'sumangle_cache');
sa_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'align', 'offset', 'delay_after', delay_after, 'cache', 'sumangle_cache');

pc_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine], 'cache', 'pixel_cache');
pc_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'align', 'offset', 'delay_after', delay_after, 'cache', 'pixel_cache');


%% 1. Setup Parameters
a = 50; 
% Updated thresholds with 70 included, sorted descending for logical fading
threshold_d = [80, 70, 60, 50, 40, 0]; 
num_thresh = numel(threshold_d);

calc_theta = @(d) rad2deg(2 * atan(a ./ (2 * d)));

datasets_str = {'pc', 'sm', 'ss', 'sa', 'distance'};
y_labels = {'Pixel Change', 'Social Motion', 'Social Speed', 'Angular Size', 'Min Dist'};
col = cmapper; 
xlims = [-420 60];

% 2. Initialize Figure
fh = figure('color','w','Position',[100 100 1200 700]); 
tlo = tiledlayout(2, 3 , 'TileSpacing', 'compact', 'Padding', 'compact');

% 3. The Loop
for i = 1:numel(datasets_str)
    ax(i) = nexttile;
    hold on
    leg_h = [];
    
    data_var_name = [datasets_str{i} '_freeze'];
    raw_data = eval(data_var_name); 
    
    % Base color for this specific tile
    base_color = col.var.(datasets_str{i});
    
    for j = 1:num_thresh
        curr_d = threshold_d(j);
        
        % 1. Masks
        [~, contact_mask] = impose_contact_threshold(bouts_proc, 'threshold', curr_d);
        total_mask = ~contact_mask & (bouts_proc.durations <= 630);
        current_data = raw_data(total_mask, :);
        
        % 2. Stats
        t = (-size(current_data, 2) + 1 : 1 : 0) + delay_after;
        avg = mean(current_data, 1, 'omitnan');
        sem = std(current_data, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(current_data), 1));
        
        % 3. Calculate Dynamic Color
        % As j increases (d decreases/changes), we move from base_color toward white
        % logic: color = base_color + (1 - base_color) * factor
        % factor 0 = full color, factor 1 = white
        if curr_d == 0
            line_col = [0 0 0]; % Keep d=0 black for contrast
            fade_factor = 0.2;
        else
            % Creates a gradient from dark/saturated to light/faded
            fade_factor = (num_thresh - j) / (num_thresh);
            line_col = base_color + (1 - base_color) * (fade_factor * 0.7);
        end
        
        % 4. Plot Shaded SEM
        fill_x = [t, fliplr(t)];
        fill_y = [avg + sem, fliplr(avg - sem)];
        fill(fill_x, fill_y, line_col, 'EdgeColor', 'none', 'FaceAlpha', 0.15, 'HandleVisibility', 'off');
        
        % 5. Plot Mean Line
        % All lines are now solid ('-'), but distinguished by color intensity
        p = plot(t, avg, 'LineWidth', 2, 'Color', line_col, ...
            'DisplayName', sprintf('d=%d px (%d%%)', curr_d, round(100 * sum(~contact_mask)/numel(contact_mask))));
                
    end
    
    % Formatting
    apply_generic(ax(i), 'font_size', 18, 'tick_length', 0.015, 'line_width', 2, ...
        'xticks', xlims(1):60:xlims(2), 'xlim', xlims);
    
    % "Offset" x-label
    current_ticks = ax(i).XAxis.TickValues;
    labels = arrayfun(@num2str, current_ticks, 'UniformOutput', false);
    labels(current_ticks == 0) = {'Offset'};
    xticklabels(ax(i), labels)
    
    ylabel(y_labels{i})
    xline(0, '--k', 'HandleVisibility', 'off');
    xline(-30, '-.k', 'HandleVisibility', 'off');
    xline(+30, '-.k', 'HandleVisibility', 'off');

    if i == 1
        legend(leg_h, 'Location', 'northwest', 'FontSize', 14, 'Interpreter', 'none');
    end
end

linkaxes(ax, 'x')

exporter(fh, paths, 'signals_multi_threshold.pdf')
exporter(fh, paths, 'signals_multi_threshold.png')


