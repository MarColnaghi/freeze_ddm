
clear all

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 3; threshold_mob = 3; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_freezes/le', 'bouts_id', id_code);
load(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.dataset, 'motion_cache.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le');
points.censoring = 10.5;
points.truncation = 0.5;

% Extract Social Motion TimeSeries
chunk_len = points.censoring * 60;

for idx_trials = 1:height(bouts_proc)

    ons = bouts_proc.onsets(idx_trials);
    sum_motion = motion_cache(bouts_proc.fly(idx_trials));
    sm_raw{idx_trials} = sum_motion(ons:ons + chunk_len - 1);

end

soc_mot_array = cell2mat(sm_raw)';

%  Now we added our vector column to the bouts table
figure
histogram(bouts_proc.avg_fs_1s_norm)
fs_quant = [0 0.45 0.9 1.35 2];
xline(fs_quant)

%%
estimates = nan(4, 2, 3);
for idx_quantiles = 1:4

    for idx_ln = 1:2
        mask = bouts_proc.sloom_norm == idx_ln & bouts_proc.avg_fs_1s_norm > fs_quant(idx_quantiles) & bouts_proc.avg_fs_1s_norm < fs_quant(idx_quantiles + 1);
        bouts_quant = bouts_proc(mask, :);
        extra.soc_mot_array = soc_mot_array(mask, :);
        model_results = run_fitting_newer_bads_only(bouts_quant, points, 'ed1', paths, 'export', false, 'extra', extra);
        estimates(idx_quantiles, idx_ln, :) = model_results.starting_position;

    end

end



%%

fh = figure('color','w','Position',[100,100, 600, 270]);
tiledlayout(1, 3, 'TileSpacing', 'loose', 'Padding', 'compact')


for idx_params = 1:size(estimates, 3)
    nexttile
    hold on
    scatter(1:4, estimates(:, 1, idx_params), 50, 'r', 'filled')
    scatter(1:4, estimates(:, 2, idx_params), 50, 'b', 'filled')
    axis square

    xlim([0 5])
    xticks([1 2 3 4])
    if idx_params == 2
        xlabel('Ranges of Focal Speed')
    end
    if idx_params == 1
        ylim([0 8])
    elseif idx_params == 2
        ylim([0 0.5])

    elseif idx_params == 3
        ylim([0 0.1])

    end
    
    apply_generic(gca)
end

exporter(fh, paths, 'quantiles.pdf')
%%
bouts_quant.sm = bouts_quant.avg_sm_freeze_norm;
bouts_quant.smp = bouts_quant.avg_sm_freeze_norm;
bouts_quant.fs = bouts_quant.avg_fs_freeze_norm;
bouts_quant.ln = bouts_quant.nloom_norm;
bouts_quant.ls = bouts_quant.sloom_norm;
bouts_quant.intercept = ones(height(bouts_quant), 1);
[nll, f, fd] = nll_fly_ddm_newer(model_results.starting_position, bouts_quant, points, 'model_ed1', 'iid', 'p', extra);

figure
histogram(bouts_quant.durations_s, 0:1/60:11, 'Normalization', 'pdf')
hold on
plot(fd, f * 60, 'k--')
%%
extra.soc_mot_array = soc_mot_array;
model_results = run_fitting_newer_bads_only(bouts_proc, points, 'ed5', paths, 'export', false, 'extra', extra);
%%
bouts_proc.sm = bouts_proc.avg_sm_freeze_norm;
bouts_proc.smp = bouts_proc.avg_sm_freeze_norm;
bouts_proc.fs = bouts_proc.avg_fs_freeze_norm;
bouts_proc.ln = bouts_proc.nloom_norm;
bouts_proc.ls = bouts_proc.sloom_norm;
bouts_proc.intercept = ones(height(bouts_proc), 1);
[nll, f, fd] = nll_fly_ddm_newer(model_results.starting_position, bouts_proc, points, 'model_ed5', 'iid', 'p', extra);

%%
fh = figure('color','w','Position',[100,100, 600, 400]);
histogram(bouts_proc.durations_s, 0:1/10:11, 'Normalization', 'pdf', 'EdgeColor', 'none')
hold on
plot(fd, f, 'k--', 'LineWidth', 2)
xlabel('Freeze Duration (s)')
ylabel('Density')
ylim([-0.02 1.52])
apply_generic(gca)
exporter(fh, paths, 'fitsss.pdf')

