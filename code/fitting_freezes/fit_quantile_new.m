
% Load the table first. We will take advantage of an already existing
% dataset.
col = cmapper();
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_freezes/le', 'bouts_id', id_code);
load(fullfile(paths.dataset, 'bouts.mat'));
thresholds = define_thresholds;
thresholds.le_window_fl = [5 40];
thresholds.le_window_sl = [15 50];
bouts = bouts_formatting(bouts, thresholds);
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le');
points.censoring = 10.5;
points.truncation = min(bouts_proc.durations_s);
link_logistic = @(x) 1./(1 + exp(-x));

motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

model_2_fit = 'ded0';

kde_estimates = importdata(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/fitting_freezes/bsl/kde_spontaneous', id_code, 'kde_estimates_bsl.mat'));
[~,idx] = unique(kde_estimates.Fkde, 'last');
extra.Fkde = kde_estimates.Fkde(idx); extra.xkde = kde_estimates.xkde(idx); extra.fkde = kde_estimates.fkde(idx); 

% Extract Social Motion TimeSeries
chunk_len = points.censoring * 60;

for idx_trials = 1:height(bouts_proc)

    ons = bouts_proc.onsets(idx_trials);
    sum_motion = motion_cache(bouts_proc.fly(idx_trials));
    sm_raw{idx_trials} = sum_motion(ons:ons + chunk_len - 1);

end

extra.soc_mot_array = cell2mat(sm_raw)';

%  Now we added our vector column to the bouts table
fh = figure('Position', [100 100 500 400], 'Color', 'w');
bouts_proc = bouts_proc(bouts_proc.durations_s >= points.truncation, :);
tiledlayout(2, 1, 'TileSpacing', 'loose')
nexttile
histogram(bouts_proc.avg_fs_1s_norm, 0:0.05:4, 'FaceColor', col.vars.fs(round(end/2), :), 'EdgeColor', 'none')
fs_quant = prctile(bouts_proc.avg_fs_1s_norm, [0, 25, 50, 75, 100]); fs_quant(1) = 0; fs_quant(end) = 2; 
fs_quant = [0 0.45 0.7 1.1 2.2]; 
xline(fs_quant);
apply_generic(gca)
nexttile
histogram(bouts_proc.avg_sm_freeze_norm, 0:0.05:4, 'FaceColor', col.vars.sm(round(end/2), :), 'EdgeColor', 'none');
sm_quant = prctile(bouts_proc.avg_sm_freeze_norm, [0, 25, 50, 75, 100]); sm_quant(1) = 0; sm_quant(end) = 2; 
sm_quant = [0, .2, 0.45, 0.8, 1.85]; 
xline(sm_quant);
apply_generic(gca)


%%
% Now you should fit a model for each quantile
% of SOCIAL MOTION

n_quantiles = 4;
n_looms = length(unique(bouts_proc.sloom_norm));
n_params = 6;
estimates = nan(n_quantiles, n_looms, n_params);
x_sm = nan(n_quantiles, n_looms);

for idx_quantiles = 1:n_quantiles

    for idx_ln = unique(bouts_proc.sloom_norm)'
        
        mask = bouts_proc.sloom_norm == idx_ln & bouts_proc.avg_sm_freeze_norm >= sm_quant(idx_quantiles) & bouts_proc.avg_sm_freeze_norm < sm_quant(idx_quantiles + 1);
        bouts_quant = bouts_proc(mask, :);
        
        figure
        histogram(bouts_quant.durations_s, [-1/120:1/10:points.censoring])
        x_sm(idx_quantiles, idx_ln) = median(bouts_quant.avg_sm_freeze_norm);

        % Here you define the Extras
        % extra = [];
        % extra.soc_mot_array = soc_mot_array(mask, :);

        model_results = run_fitting_newer(bouts_quant, points, model_2_fit, paths, 'export', false, 'extra', extra);
        est = table2array(model_results.estimates_mean);

        estimates(idx_quantiles, idx_ln, :) = est(~isnan(est));

    end
end

estimates( :, :, 5) = link_logistic(estimates( :, :, 5));

fh = figure('color','w','Position',[100,100, 700, 320]);
tiledlayout(1, n_params, 'TileSpacing', 'loose', 'Padding', 'compact')

for idx_params = 1:size(estimates, 3)
    nexttile
    hold on
    scatter(mean(x_sm, 2), estimates(:, 1, idx_params), 50, 'r', 'filled')
    scatter(mean(x_sm, 2), estimates(:, 2, idx_params), 50, 'b', 'filled')
    axis square

    xlim([0 1.5])
    xticks(mean(x_sm, 2))
    xticklabels([1 2 3 4])

    if idx_params == 1
        ylim([0 2])

    elseif idx_params == 2
        ylim([0 1.5])

    elseif idx_params == 3
        ylim([0 1])

    elseif idx_params == 4
        ylim([0 0.2])

    end
    
    apply_generic(gca)
end

%% Now you should fit a model for each quantile
% of FOCAL SPEED

n_quantiles = 4;
n_looms = length(unique(bouts_proc.sloom_norm));
n_params = 6;
estimates = nan(n_quantiles, n_looms, n_params);
x_fs = nan(n_quantiles, n_looms);

for idx_quantiles = 1:n_quantiles

    for idx_ln = unique(bouts_proc.sloom_norm)'

        mask = bouts_proc.sloom_norm == idx_ln & bouts_proc.avg_fs_1s_norm >= fs_quant(idx_quantiles) & bouts_proc.avg_fs_1s_norm < fs_quant(idx_quantiles + 1);
        ec.soc_mot_array = extra.soc_mot_array(mask, :);
        
        bouts_quant = bouts_proc(mask, :);

        x_fs(idx_quantiles, idx_ln) = median(bouts_quant.avg_fs_1s_norm);

        % Here you define the Extras
        % extra = [];
        % extra.soc_mot_array = soc_mot_array(mask, :);

        model_results = run_fitting_newer(bouts_quant, points, model_2_fit, paths, 'export', false, 'extra', ec);
        est = table2array(model_results.estimates_mean);

        estimates(idx_quantiles, idx_ln, :) = est(~isnan(est));

    end
end

estimates(:, :, 5) = link_logistic(estimates( :, :, 5));
estimates(3, :, [1 2 3 4]) = estimates(3, :, [3 4 1 2])
estimates(3, :, 5) = 1 - estimates(3, :, 5)

fh = figure('color','w','Position',[100,100, 700, 320]);
tiledlayout(1, n_params, 'TileSpacing', 'loose', 'Padding', 'compact')


for idx_params = 1:size(estimates, 3)
    nexttile
    hold on
    scatter(mean(x_fs, 2), estimates(:, 1, idx_params), 50, 'r', 'filled')
    scatter(mean(x_fs, 2), estimates(:, 2, idx_params), 50, 'b', 'filled')
    axis square

    xlim([0 2])
    xticks(mean(x_fs, 2))
    xticklabels([1 2 3 4])

    if idx_params == 1
        ylim([0 2])
    elseif idx_params == 2
        ylim([0 1.5])

    elseif idx_params == 3
        ylim([0 1])

    elseif idx_params == 4
        ylim([0 0.2])

    elseif idx_params == 5
        ylim([0 1])
    end
    
    apply_generic(gca)
end


% Now you should fit a model for each quantile
% of LOOM NUMBER
