dx = 1/60;
xbin = (0:dx:dx*(ceil(10/dx)+10))';
xxi     = 0:dx:1200;
h       = .1;

threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', fullfile('fitting_freezes','le'), 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));

thresholds = define_thresholds;
% thresholds.le_window_fl = [5 40];
% thresholds.le_window_sl = [15 50];

bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);
freezes = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);

fh = figure('Position', [100 100 500 600], 'Color', 'w');
t = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
i = 0;

for idx_sm = 1
    for idx_ls = 2
        for idx_fs = 1

            i = i + 1;
            nexttile(t)
            hold on

            ax(i) = gca;

            [freezes_quant, ~] = quantilizer(freezes, 'idx_quanti', struct('sm', idx_sm, 'ls', idx_ls, 'fs', idx_fs), 'quant', struct('sm', 1, 'ls', 2, 'fs', 2));

            RTs{1,1} = freezes_quant.durations_s;
            RTs{2,1} = freezes_quant.durations_s;
            % RTD = kreg_single(RTs, RTs, xxi, xbin, h, 0, 500);

            censored = freezes_quant.durations_s > 10.5;
            [ks_noc] = ksdensity(freezes_quant.durations_s, xxi, 'BoundaryCorrection', 'reflection', 'Bandwidth', h, 'Support', [min(freezes.durations_s) - 0.001, max(freezes.durations_s) + 0.001 ]);
            [ks] = ksdensity(freezes_quant.durations_s, xxi, 'BoundaryCorrection', 'reflection', 'Bandwidth', h, 'Support', [min(freezes.durations_s) - 0.001, max(freezes.durations_s) + 0.001 ], 'Censoring', censored);

            plot(xxi, ks_noc, 'LineWidth', 2, 'Color', [0.8 0.8 0.8])
            histogram(freezes_quant.durations_s, min(freezes_quant.durations_s):10/60:600, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8])
        end
    end
end

apply_generic(gca, 'xlim', [0 10.5])
xlabel('Duration (s)')