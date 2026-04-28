clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion/clustering', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
%thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);

%% Extract the pre-freezing social motion
sm_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat');
inizio = 0;
fine = 630;
fps = 60;
t = (inizio:fine) / fps;
sm_freeze_full = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine]);

%%
fh = figure('Position', [100 100 1300 700], 'Color', 'w');
tl = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
i = 0;
bin_size = 0.1;
for idx_sm = 1:3
    for idx_ls = 0:1
        for idx_fs = 1:2

            nexttile
            hold on
            i = i + 1;

            [freezes_quant, mask] = quantilizer_v2(bouts_proc, 'indexed_quantile', struct('sm', idx_sm, 'ls', idx_ls, 'fs', idx_fs));

            sum_freezes(i) = numel(mask);
            histogram(freezes_quant.durations_s, min(freezes_quant.durations_s):bin_size:max(freezes_quant.durations_s), 'Normalization', 'pdf')
            apply_generic(gca, 'ylim', [0 2], 'xlim', [0 10.5])

            data = sm_freeze(mask, :);

            mu  = mean(data, 1, 'omitnan');
            sem = std(data, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data), 1));

            t_full = (1:numel(mu)) ./ fps;
            m = ~isnan(sem);
            
            yyaxis right
            % --- SHADED SEM ---
            fill_between(t_full(m), ...
                mu(m) + sem(m), ...
                mu(m) - sem(m), ...
                [], ...
                'FaceColor', [0.8 0.24 0.24], ... % or use a color per condition
                'EdgeColor', 'none', ...
                'FaceAlpha', 0.2);

            % --- MEAN TRACE ---
            plot(t_full, mu, 'Color', [0.8 0.24 0.24], 'LineWidth', 1.5)

            ylim([0 1])

        end
    end
end