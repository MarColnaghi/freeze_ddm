clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'collision_analysis', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [0 55], 'le_window_fl', [0 55]));
% thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);

thresholds_collisions = 30:10:100;
tot_n_collisions_of = nan(numel(thresholds_collisions), 1);
tot_n_collisions_ili = nan(numel(thresholds_collisions), 1);

for idx_collision_threshold = 1:numel(thresholds_collisions)
    
    curr_thre = thresholds_collisions(idx_collision_threshold);
    [bouts_proc, contact_mask_of] = impose_contact_threshold(bouts_proc, 'threshold', curr_thre, 'type', 'onlyfreeze');
    tot_n_collisions_of(idx_collision_threshold) = sum(contact_mask_of);

    [bouts_proc, contact_mask_ili] = impose_contact_threshold(bouts_proc, 'threshold', curr_thre, 'type', 'ili');
    tot_n_collisions_ili(idx_collision_threshold) = sum(contact_mask_ili);

end

fh = figure('color','w','Position', [100, 100, 700, 600]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact')

nexttile
bar(thresholds_collisions, tot_n_collisions_of)
apply_generic(gca, 'font_size', 24, 'tick_length', 0.015, 'line_width', 2, 'ylim', [0 height(bouts_proc)], 'no_x', true);

xtickangle(0)
xlabel('Threshold Value')
ylabel('Count')

text(thresholds_collisions, ...
     tot_n_collisions_of + 40, ...
     string(round(tot_n_collisions_of/ height(bouts_proc), 2)), ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','bottom', 'FontSize', 16)

nexttile
bar(thresholds_collisions, tot_n_collisions_ili)
apply_generic(gca, 'font_size', 24, 'tick_length', 0.015, 'line_width', 2, 'ylim', [0 height(bouts_proc)]);
xtickangle(0)
xlabel('Threshold Value')
ylabel('Count')

text(thresholds_collisions, ...
     tot_n_collisions_ili + 40, ...
     string(round(tot_n_collisions_ili/ height(bouts_proc), 2)), ...
     'HorizontalAlignment','center', ...
     'VerticalAlignment','bottom', 'FontSize', 16)

exporter(fh, paths, 'number_of_collisions.pdf')
