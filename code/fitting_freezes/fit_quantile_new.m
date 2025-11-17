
% Load the table first. We will take advantage of an already existing
% dataset.
col = cmapper();
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_freezes/le/quantiles', 'bouts_id', id_code);
load(fullfile(paths.dataset, 'bouts.mat'));
thresholds = define_thresholds;
thresholds.le_window_fl = [1 60];
thresholds.le_window_sl = [1 60];
bouts = bouts_formatting(bouts, thresholds);
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le');
points.censoring = 10.5;
points.truncation = 0;
link_logistic = @(x) 1./(1 + exp(-x));

motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

model_2_fit = 'simed0';

kde_estimates = importdata(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/fitting_freezes/bsl/kde_spontaneous', id_code, 'kde_estimates_bsl.mat'));
[~,idx] = unique(kde_estimates.Fkde, 'last');
extra.Fkde = kde_estimates.Fkde(idx); extra.xkde = kde_estimates.xkde(idx); extra.fkde = kde_estimates.fkde(idx); 

% Extract Social Motion TimeSeries
chunk_len = points.censoring * 60;

for idx_trials = 1:height(bouts_proc)

    ons = bouts_proc.onsets(idx_trials);
    sum_motion = motion_cache(bouts_proc.fly(idx_trials)) ./ 10;
    sm_raw{idx_trials} = sum_motion(ons:ons + chunk_len - 1);

end

extra.soc_mot_array = cell2mat(sm_raw)';

%  Now we added our vector column to the bouts table
fh = figure('Position', [100 100 500 400], 'Color', 'w');
bouts_proc = bouts_proc(bouts_proc.durations_s >= points.truncation, :);
tiledlayout(2, 1, 'TileSpacing', 'loose')
nexttile
histogram(bouts_proc.avg_fs_1s_norm, 0:0.05:4, 'FaceColor', col.vars.fs(round(end/2), :), 'EdgeColor', 'none')
quant.avg_fs_1s_norm = prctile(bouts_proc.avg_fs_1s_norm, [0, 25, 50, 75, 100]); quant.avg_fs_1s_norm(1) = 0; quant.avg_fs_1s_norm(end) = 2; 
quant.avg_fs_1s_norm = [0 0.45 0.7 1.0 2.2]; 
xline(quant.avg_fs_1s_norm);
apply_generic(gca)
nexttile
histogram(bouts_proc.avg_sm_freeze_norm, 0:0.05:4, 'FaceColor', col.vars.sm(round(end/2), :), 'EdgeColor', 'none');
quant.avg_sm_freeze_norm = prctile(bouts_proc.avg_sm_freeze_norm, [0, 25, 50, 75, 100]); quant.avg_sm_freeze_norm(1) = 0; quant.avg_sm_freeze_norm(end) = 2; 
quant.avg_sm_freeze_norm = [0, .2, 0.45, 0.8, 1.85]; 
xline(quant.avg_sm_freeze_norm);
apply_generic(gca)

quant.nloom_norm = [-0.05 0.55 1.05 1.55 2.05];

%% Now you should fit a model for each quantile
% of SOCIAL MOTION

n_quantiles = 4;
n_looms = length(unique(bouts_proc.sloom_norm));
n_params = 2;
estimates = nan(n_quantiles, n_looms, n_params);
x = nan(n_quantiles, n_looms);
vars = {'avg_sm_freeze_norm', 'avg_fs_1s_norm', 'nloom_norm'};

paths = path_generator('folder', fullfile('fitting_freezes/le/quantiles', model_2_fit), 'bouts_id', id_code);
create_output_dirs(paths)

for idx_vars = vars
    vr = idx_vars{1};

    for idx_quantiles = 1:n_quantiles

        for idx_ln = unique(bouts_proc.sloom_norm)'

            edges = quant.(vr);
            mask = bouts_proc.sloom_norm == idx_ln & bouts_proc.(vr) >= edges(idx_quantiles) & bouts_proc.(vr) < edges(idx_quantiles + 1);
            ec.soc_mot_array = extra.soc_mot_array(mask, :);


            bouts_quant = bouts_proc(mask, :);
            x(idx_quantiles, idx_ln) = median(bouts_quant.(vr));

            figure
            histogram(bouts_quant.durations_s, -1/120:1/5:10.5, 'Normalization', 'pdf')
            drawnow

            model_results = run_fitting_newer(bouts_quant, points, model_2_fit, paths, 'export', false, 'extra', ec, 'pass_ndt', true);
            model_results.quant = quant.(vr);

            est = table2array(model_results.estimates_mean);

            estimates(idx_quantiles, idx_ln, :) = est(~isnan(est));

        end
    end

    % estimates( :, :, 5) = link_logistic(estimates( :, :, 5));
    mkdir(fullfile(paths.results, (vr))); cd(fullfile(paths.results, (vr)));

    save('estimates.mat', 'estimates')
    save('model_results.mat', 'model_results');


    % Now you should fit a model for each quantile
    % of FOCAL SPEED


    % estimates(:, :, 5) = link_logistic(estimates( :, :, 5));
    % estimates(3, :, [1 2 3 4]) = estimates(3, :, [3 4 1 2])
    % estimates(3, :, 5) = 1 - estimates(3, :, 5)

    fh = figure('color','w','Position',[100,100, 700, 320]);
    tiledlayout(1, n_params, 'TileSpacing', 'loose', 'Padding', 'compact')


    for idx_params = 1:size(estimates, 3)
        nexttile
        hold on
        scatter(mean(x, 2), estimates(:, 1, idx_params), 50, 'r', 'filled')
        scatter(mean(x, 2), estimates(:, 2, idx_params), 50, 'b', 'filled')
        axis square

        xlim([0 2])
        xticks(mean(x, 2))
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
end


% Now you should fit a model for each quantile
% of LOOM NUMBER
function create_output_dirs(paths)
    % Ensure base directories exist
    if ~exist(paths.fig, 'dir'), mkdir(paths.fig); end
    if ~exist(paths.results, 'dir'), mkdir(paths.results); end

    % Auto-incrementing run folder inside results
    run_folders = dir(fullfile(paths.results, 'run*'));
    run_nums = [];

    for i = 1:length(run_folders)
        if run_folders(i).isdir
            tokens = regexp(run_folders(i).name, '^run(\d+)$', 'tokens');
            if ~isempty(tokens)
                run_nums(end+1) = str2double(tokens{1}{1}); %#ok<AGROW>
            end
        end
    end

    if isempty(run_nums)
        next_run = 1;
    else
        next_run = max(run_nums) + 1;
    end

    run_name = sprintf('run%02d', next_run);
    paths.results = fullfile(paths.results, run_name);
    mkdir(paths.results);

    % Also update figure path to match the new run
    paths.fig = fullfile(paths.fig, run_name);
    mkdir(paths.fig);

    % Assign the updated paths back to base workspace (if needed)
    assignin('caller', 'paths', paths);
end