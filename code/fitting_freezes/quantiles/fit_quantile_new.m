clearvars

col = cmapper();
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_freezes/le/quantiles', 'bouts_id', id_code);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

thresholds = define_thresholds;
%thresholds.le_window_fl = [1 60];
%thresholds.le_window_sl = [1 60];
% thresholds.le_window_fl = [0 45];
% thresholds.le_window_sl = [10 55];

bouts = bouts_formatting(bouts, thresholds);
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);

total_length = 30;

[mean_sm_before_freeze, mean_sm_during_freeze] = extract_sm_columns(bouts_proc, motion_cache, 'chunk_dur', total_length);
bouts_proc.smp = mean_sm_before_freeze;

% At some point we should modify this piece of code and unify
% extract_sm_columns and extract_sm_from_bouts.

% sm_before_freeze = mean(extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'window', [-31 -1]), 2);

points.censoring = 10.5;
points.truncation = 0.5;

%extra.soc_mot_array = cell2mat(sm_during)';

%  Now we added our vector column to the bouts table
bouts_proc = bouts_proc(bouts_proc.durations_s >= points.truncation, :);

fh = figure('Position', [100 100 500 400], 'Color', 'w');
tiledlayout(2, 1, 'TileSpacing', 'loose')

nexttile
histogram(bouts_proc.fs, 0:0.02:3, 'FaceColor', col.vars.fs(round(end/2), :), 'EdgeColor', 'none')
quant.fs = prctile(bouts_proc.fs, [0, 25, 50, 75, 100]); quant.fs(1) = 0; quant.fs(end) = 2;
quant.fs = [0 0.45 0.7 1.0 2.2];
xline(quant.fs);
apply_generic(gca, 'xlim', [0 2])

nexttile
histogram(bouts_proc.sm, 0:0.02:3, 'FaceColor', col.vars.sm(round(end/2), :), 'EdgeColor', 'none');
quant.sm = prctile(bouts_proc.sm, [0, 25, 50, 75, 100]); quant.sm(1) = 0; quant.sm(end) = 1.2;
quant.sm = [0, .25, 0.45, 0.8, 1.85];
xline(quant.sm);
apply_generic(gca, 'xlim', [0 2])

quant.nloom_norm = [-0.05 0.55 1.05 1.55 2.05];

%% Now you should fit a model for each quantile
% of SOCIAL MOTION

model_2_fit = 'dddim0';
model = str2func(strcat('model_', model_2_fit));
model_struct = model();

n_quantiles = 3;
p = linspace(0, 100, n_quantiles + 1);

quant.sm = prctile(bouts_proc.sm, p);
n_looms = length(unique(bouts_proc.ls));
n_params = length(fieldnames(model_struct));

estimates = nan(n_quantiles, n_looms, n_params);
x = nan(n_quantiles, n_looms);
% vars = {'sm', 'fs', 'nloom_norm'};
vars = {'sm'};

paths = path_generator('folder', fullfile('fitting_freezes/le/quantiles', model_2_fit), 'bouts_id', id_code);
create_output_dirs(paths)
i = 0;

for idx_vars = vars
    vr = idx_vars{1};

    for idx_quantiles = 1:n_quantiles

        for idx_ln = unique(bouts_proc.ls)'

            i = i + 1;

            edges = quant.(vr);
            mask = bouts_proc.ls == idx_ln & bouts_proc.(vr) >= edges(idx_quantiles) & bouts_proc.(vr) < edges(idx_quantiles + 1);


            bouts_quant = bouts_proc(mask, :);
            x(idx_quantiles, idx_ln) = median(bouts_quant.(vr));

            fh(i) = figure;
            histogram(bouts_quant.durations_s, min(bouts_quant.durations_s):1/20:10.5, 'Normalization', 'pdf')
            hold on

            [~ , estimate] = run_fitting_newer(bouts_quant, points, model_2_fit, paths, 'export', false, 'only_bads', true);

            estimates(idx_quantiles, idx_ln, :) = estimate;

            bouts_quant.intercept = ones(height(bouts_quant), 1);

            [~, f, fd] = nll_fly_ddm_newer(estimate, bouts_quant, points, func2str(model), 'iid','p', []);

            figure(fh(i))
            plot(fd, f, 'r--')

            drawnow

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
            ylim([0 3.5])
        elseif idx_params == 2
            ylim([0 1.5])

        elseif idx_params == 3
            ylim([0 2])

        elseif idx_params == 4
            ylim([0 5])

        elseif idx_params == 5
            ylim([0 10])
        end

        apply_generic(gca)
    end
end

%%
fh = figure('color','w','Position',[100,100, 550, 350]);
tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'compact')


for idx_params = 1:2
    nexttile
    hold on
    xvals = mean(x, 2);          % should be 4x1
    yvals = estimates(:, 2, idx_params);

    % build color matrix (one row per point)
    C = repmat([0.7 0.7 0.7], numel(xvals), 1);

    for q = 1:4
        C(q, :) = col.vars.sm(1 + q, :);
    end

    scatter(mean(x, 2), estimates(:, 1, idx_params), 140, 'r', 'filled')
    scatter(mean(x, 2), estimates(:, 2, idx_params), 140, 'k', 'filled')
    axis square

    xlim([0 2])
    xticks(mean(x, 2))
    xticklabels({'Q1', 'Q2', 'Q3', 'Q4'})

    if idx_params == 1
        ylim([-0. 1])
        yticks([0 0.5 1])
        ylabel('Evidence Strength ($\mu$)', 'Interpreter', 'latex', 'Color', col.param.mu)

    elseif idx_params == 2
        ylim([0 2.0])
        ylabel('Bound Distance ($\theta$)', 'Interpreter', 'latex', 'Color', col.param.theta)

    end

    apply_generic(gca, 'xlim', [0 1])

end

exporter(fh, paths, 'quartile_fitting_with_ls.pdf')




%%

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