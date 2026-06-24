clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion/proximity', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));
speed_cache = importdata(fullfile(paths.cache_path, 'speed_cache.mat'));
pixel_cache = importdata(fullfile(paths.cache_path, 'pixel_cache.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [0 55], 'le_window_fl', [0 55]));
% thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);
inizio = 0;
fine = 630;

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);
bouts_proc = bouts_proc(bouts_proc.durations <= fine, :);

distance_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine], 'cache', 'mindist_cache');
distance_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'align', 'onset', 'delay_after', 0, 'cache', 'mindist_cache');

threshold = 120; 
[bouts_proc, contact_mask] = impose_contact_threshold(bouts_proc, 'threshold', threshold, 'type', 'onlyfreeze');

matrix = distance_freeze(:, (inizio + 1:fine));
[value_min_mindistance, when_min_mindistance] = min(matrix, [], 2, 'omitnan');
tot_duration = (sum(~isnan(matrix), 2));

relative_when = when_min_mindistance./tot_duration;
figure
scatter(relative_when, value_min_mindistance)

hold on

scatter(relative_when(~contact_mask), value_min_mindistance(~contact_mask))
mask = relative_when > 0.2 & relative_when < 0.8;

bouts_with_contacts_during_freeze = bouts_proc(contact_mask & mask, :);
writetable(bouts_with_contacts_during_freeze, sprintf('bouts_with_contacts_during_freeze_%d.csv', threshold))

figure
histogram(value_min_mindistance(relative_when > 0.95), 0:3:600, 'Normalization', 'pdf')
hold on
histogram(value_min_mindistance(relative_when < 0.05), 0:3:600, 'Normalization', 'pdf')

distance_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'align', 'offset', 'delay_after', 0, 'cache', 'mindist_cache');
figure
histogram(mean(distance_freeze(:, end:end), 2, 'omitnan'), 0:3:600, 'Normalization', 'pdf')
xline(70, 'k-',  'LineWidth', 2)


break_with_contact = mean(distance_freeze(:, end:end), 2, 'omitnan') < threshold;
break_without_contact = mean(distance_freeze(:, end:end), 2, 'omitnan') >= threshold;

size_vector = 180;
speed_matrix = nan(height(bouts_proc), 2 * size_vector);
pixel_matrix = nan(height(bouts_proc), 2 * size_vector);

for idx_bout = 1:height(bouts_proc)

    s = bouts_proc.ends(idx_bout) - size_vector;
    e = bouts_proc.ends(idx_bout) + size_vector - 1;

    fly = bouts_proc.fly(idx_bout);

    fly_speed = speed_cache(fly);
    fly_pixel = pixel_cache(fly);

    % Make sure indices stay within bounds
    valid_e = min(e, length(fly_speed));

    % Number of available samples
    n_available = valid_e - s + 1;

    % Only fill available portion
    if s <= length(fly_speed) && n_available > 0
        speed_matrix(idx_bout, 1:n_available) = fly_speed(s:valid_e);
        pixel_matrix(idx_bout, 1:n_available) = fly_pixel(s:valid_e);

    end


end

next_bout = bouts_proc.id + 1;

% Find corresponding row in bouts
[found, idx] = ismember(next_bout, bouts.id);

% Initialize with NaN
next_duration = nan(height(bouts_proc),1);

% Keep only valid matches from the same fly
valid = found & ...
        (bouts.fly(idx) == bouts_proc.fly);

% Extract durations only for valid next bouts
next_duration(valid) = bouts.durations(idx(valid));

figure
hold on
plot(-size_vector:size_vector - 1, mean(speed_matrix, 1, 'omitnan'), 'k-')

plot(-size_vector:size_vector - 1, mean(speed_matrix(break_with_contact,:), 1, 'omitnan'), 'r-')
plot(-size_vector:size_vector - 1, mean(speed_matrix(break_without_contact,:), 1, 'omitnan'), 'g-')


fh = figure('color','w','Position', [100, 100, 700, 600]);
tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact')
nexttile(1)
ax(1) = gca;
hold on
histogram(next_duration(break_without_contact), 0:1:500, 'FaceColor', 'g', 'Normalization', 'pdf', 'EdgeColor', 'none')
histogram(next_duration(break_with_contact), 0:1:500, 'FaceColor', 'r', 'Normalization', 'pdf', 'EdgeColor', 'none')
apply_generic(gca, 'font_size', 20, 'tick_length', 0.015, 'line_width', 2);
xlabel('Duration next moving bout (frames)')

nexttile(2, [2,1])
ax(2) = gca;
hold on
plot(-size_vector:size_vector - 1, mean(pixel_matrix, 1, 'omitnan'), 'k-', 'LineWidth', 2)
plot(-size_vector:size_vector - 1, mean(pixel_matrix(break_with_contact,:), 1, 'omitnan'), 'r-', 'LineWidth', 2)
plot(-size_vector:size_vector - 1, mean(pixel_matrix(break_without_contact,:), 1, 'omitnan'), 'g-', 'LineWidth', 2)
apply_generic(ax(2), 'font_size', 20, 'tick_length', 0.015, 'line_width', 2);
xlabel('Time (frames)')
ylabel('Mean Pixel Change')
linkaxes(ax(:), 'x')
xlim([-120, 120])

%%

[bouts_proc, contact_mask] = impose_contact_threshold(bouts_proc, 'threshold', 70, 'type', 'ili');
