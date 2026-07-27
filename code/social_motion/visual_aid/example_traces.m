clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion/visual_aid', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [0 55], 'le_window_fl', [0 55]));
% thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);
bouts_proc = sortrows(bouts_proc, 'sm');
threshold = 0;
[bouts_proc, ~] = impose_contact_threshold(bouts_proc, 'threshold', threshold);

inizio = 0;
fine = 630;

align = 'offset';

sm_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'norm_factor', 10, 'cache', 'motion_cache');
sm_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine], 'norm_factor', 10, 'cache', 'motion_cache');

n_trials = 50;
trials = sort(randi(height(bouts_proc), n_trials, 1));

fh = figure('color','w','Position', [100, 100, 500, 900]);
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact')

nexttile
hold on

for idx_counter = 1:n_trials
    current_bout = trials(idx_counter);

    t_ili = inizio:fine;
    plot(t_ili, sm_ili(current_bout, :) + 2 * idx_counter, 'Color', [0.8 0.8 0.8], 'LineWidth', 1)

    idx_final = numel(sm_freeze(current_bout, :));
    t_freeze = inizio:idx_final-1;
    plot(t_freeze, sm_freeze(current_bout, :) + 2 * idx_counter, 'Color', [0.2 0.2 0.2], 'LineWidth', 1.25)

end

apply_generic(gca, 'font_size', 24, 'tick_length', 0.015, 'line_width', 2, 'no_y', true, 'xlim', [inizio fine], 'xticks', 0:120:600, 'ylim', [-2 max(sm_freeze(current_bout, :) + 2 * idx_counter) + 2]);
xticklabels({'Onset', '2', '4', '6', '8', '10'});

% Set the interpreter to TeX so it recognizes the \newline command
ax = gca;
ax.TickLabelInterpreter = 'tex';
xtickangle(0);
xlabel('Time (s)')
exporter(fh, paths, 'traces_examples.pdf')
exporter(fh, paths, 'traces_examples.png')

%%
% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);
bouts_proc = sortrows(bouts_proc, 'sm');
threshold = 0;
[bouts_proc, contact_mask] = impose_contact_threshold(bouts_proc, 'threshold', threshold);

inizio = 0;
fine = 630;

align = 'offset';

sm_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'norm_factor', 10, 'cache', 'motion_cache', 'align', 'offset');
%sm_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine], 'norm_factor', 10, 'cache', 'motion_cache');

%%

n_trials = 50;
trials = sort(randi(height(bouts_proc), n_trials, 1));

t = inizio:fine;

fh = figure('color','w','Position', [100, 100, 550, 900]);
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact')

nexttile
hold on

for idx_counter = 1:n_trials
    current_bout = trials(idx_counter);

    t_freeze = -size(sm_freeze, 2):1:-1;
    plot(t_freeze, sm_freeze(current_bout, :) + 2 * idx_counter, 'Color', [0.2 0.2 0.2], 'LineWidth', 1.25)

end

apply_generic(gca, 'font_size', 24, 'tick_length', 0.015, 'line_width', 2, 'no_y', true, 'xlim', [-630 0], 'xticks', -600:120:0, 'ylim', [-2 max(sm_freeze(current_bout, :) + 2 * idx_counter) + 2]);
xticklabels({'-10', '-8', '-6', '-4', '-2', 'Offset'});

% Set the interpreter to TeX so it recognizes the \newline command
ax = gca;
ax.TickLabelInterpreter = 'tex';
xtickangle(0);
xlabel('Time (s)')
exporter(fh, paths, 'traces_examples_offset.pdf')
exporter(fh, paths, 'traces_examples_offset.png')
