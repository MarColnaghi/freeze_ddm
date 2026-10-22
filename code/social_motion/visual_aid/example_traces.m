clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 3; threshold_mob = 3; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion/visual_aid', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);
bouts_proc = sortrows(bouts_proc, 'sm');
threshold = 70;
[bouts_proc, ~] = impose_contact_threshold(bouts_proc, 'threshold', threshold);

% If set, keep only the freezes that ended on their own within the window: no
% contact, not cut at the cap, not running into the end of the recording.
remove_censored = true;
if remove_censored
    bouts_proc = bouts_proc(~bouts_proc.is_censored, :);
end

inizio = 0;
fine = 630;

sm_freeze_onset = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'norm_factor', 10, 'cache', 'motion_cache');
sm_freeze_offset = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat', 'norm_factor', 10, 'cache', 'motion_cache', 'align', 'offset');
sm_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine], 'norm_factor', 10, 'cache', 'motion_cache');

% The same bouts are shown in the onset- and offset-aligned figures. Rows
% are sorted by sm, so the stacked traces go from low to high sm.
rng(1);
n_trials = 50;
trials = sort(randperm(height(bouts_proc), n_trials));
spacing = 0.7;

%% Onset-aligned
fh = figure('color','w','Position', [100, 100, 500, 900]);
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact')

nexttile
hold on

t_ili = inizio:fine;
t_freeze = inizio:size(sm_freeze_onset, 2) - 1;
for idx_counter = 1:n_trials
    current_bout = trials(idx_counter);

    plot(t_ili, sm_ili(current_bout, :) + spacing * idx_counter, 'Color', [0.8 0.8 0.8], 'LineWidth', 1)
    plot(t_freeze, sm_freeze_onset(current_bout, :) + spacing * idx_counter, 'Color', [0.2 0.2 0.2], 'LineWidth', 1.25)
end

y_top = spacing * n_trials + max([sm_ili(trials, :), sm_freeze_onset(trials, 1:min(end, fine + 1))], [], 'all', 'omitnan') + spacing;
apply_generic(gca, 'font_size', 24, 'tick_length', 0.015, 'line_width', 2, 'no_y', true, 'xlim', [inizio fine], 'xticks', 0:120:600, 'ylim', [-spacing y_top]);
xticklabels({'Onset', '2', '4', '6', '8', '10'});
xtickangle(0);
xlabel('Time (s)')
exporter(fh, paths, 'traces_examples.pdf')
exporter(fh, paths, 'traces_examples.png')

%% Offset-aligned
fh = figure('color','w','Position', [100, 100, 450, 570]);
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact')

nexttile
hold on

t_freeze = -size(sm_freeze_offset, 2):-1;
for idx_counter = 1:n_trials
    current_bout = trials(idx_counter);

    plot(t_freeze, sm_freeze_offset(current_bout, :) + spacing * idx_counter, 'Color', [0.2 0.2 0.2], 'LineWidth', 1.25, 'Clipping', 'off')
end

y_top = spacing * n_trials + max(sm_freeze_offset(trials, max(1, end - fine):end), [], 'all', 'omitnan') + spacing;
apply_generic(gca, 'font_size', 24, 'tick_length', 0.015,c 'line_width', 2, 'no_y', true, 'xlim', [-fine 0], 'xticks', -600:120:0, 'ylim', [-spacing y_top]);
xticklabels({'-10', '-8', '-6', '-4', '-2', 'Break'});
xtickangle(0);
xlabel('Time (s)')
xlim([-400 0])
exporter(fh, paths, 'traces_examples_offset.pdf')
exporter(fh, paths, 'traces_examples_offset.png')
