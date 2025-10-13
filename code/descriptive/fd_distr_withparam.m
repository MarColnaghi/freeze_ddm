%% Script Name:
% fd_distr_with_param_new
% last update: 18.02

%% Description:
% Plots the distribution of freeze duration as a function of a specific
% experimental parameter

%% Load Data and Preprocess
%clear all; %close all; clc;

function fd_distr_withparam(bouts_proc, moment, paths, export)

% Load Colors
extra.quantiles = 5;
col = cmapper([], extra);

if nargin == 0

    disp('no input provided, will load the dataset in the dataset folder')

    paths = path_generator('descr/fd_distr');

    thresholds = define_thresholds;

    % Load the Data
    load(fullfile(paths.dataset, 'bouts.mat'))
    bouts = bouts(bouts.genotype == 1, :);
    bouts = bouts_formatting(bouts, thresholds);
    [bouts_proc] = data_parser_new(bouts);
    moment = 'loom';

    if strcmp(moment, 'bsl')
        bouts_proc = bouts_proc(bouts_proc.period == 0, :);

    elseif strcmp(moment, 'loom')
        bouts_proc = bouts_proc(bouts_proc.period == 1, :);

    end
else
    disp('input provided')
    paths.fig = fullfile(paths.fig, 'fd_durs');

end

%%%%%%%%%%%%%%%%%%%%%
% Choose Param here %
%%%%%%%%%%%%%%%%%%%%%

%param_list = {'onsets_loomaligned'}; 
param_list =  {'avg_sm_freeze_norm', 'avg_fs_1s_norm', 'sloom', 'nloom', 'moving_flies', 'onsets_loomaligned_norm'};

for idx_param = param_list

    param = idx_param{1};

    if strcmp(param, 'nloom')
        num_quantiles = 4;
        param_str = 'Loom Number';
        if strcmp(moment, 'bsl')
            thresholds = 1:5:21;
        elseif  strcmp(moment, 'loom')
            thresholds = 1:5:21;
        end

    elseif strcmp(param, 'sloom')
        num_quantiles = 2;
        param_str = 'Loom Speed';
        thresholds = 12.5:25:75;

    elseif strcmp(param, 'avg_fs_1s_norm')
        num_quantiles = 4;
        param_str = 'Focal Speed';
        thresholds = quantile(bouts_proc.(param), linspace(0, 1, num_quantiles + 1));

    elseif strcmp(param, 'avg_sm_freeze_norm')
        num_quantiles = 4;
        param_str = 'Avg. Social Motion';
        thresholds = quantile(bouts_proc.(param), linspace(0, 1, num_quantiles + 1));

    elseif strcmp(param, 'moving_flies')
        num_quantiles = 5;
        param_str = 'N. Moving Flies';
        thresholds = -0.5:1:5.5;
    
    elseif strcmp(param, 'onsets_loomaligned_norm')
        num_quantiles = 4;
        param_str = 'Onset';
        thresholds = quantile(bouts_proc.(param), linspace(0, 1, num_quantiles + 1));
    end

    quantiles = discretize(bouts_proc.(param), thresholds);

    %%%%%%%%%%%%%%%%%%%%%
    % Choose Param here %
    %%%%%%%%%%%%%%%%%%%%%

    type = 'ecdf';
    bw = 0.25;

    dx = .005;
    xbin = (0:dx:dx*(ceil(10/dx)+10))';
    x     = (0:10:10000)/1000;
    h       = .05;

    % Load Colors
    extra.quantiles = num_quantiles;
    col = cmapper([], extra);

    % Create Figure
    fh = figure('color','w','Position',[700, 100, 400, 550]);
    tl = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'loose');

    nexttile
    ax = gca;
    hold on

    for idx_quant = 1:num_quantiles
        if strcmp(type, 'ecdf')
            [f, x] = ecdf(bouts_proc.durations_s(quantiles == idx_quant), 'censoring', bouts_proc.bout_with_loom(quantiles == idx_quant));
            length(find(bouts_proc.durations_s(quantiles == idx_quant)))
            plot(x, f, 'Color', col.vars.(param)(1 + idx_quant,:), 'LineWidth', 3)
            ylabel('ecdf');
            ylim([0 1]);

        elseif strcmp(type, 'kde')
            RTs{1,1} = bouts_proc.durations_s(quantiles == idx_quant); RTs{2,1} = bouts_proc.durations_s(quantiles == idx_quant);
            [f] = kreg_single(RTs, RTs, x, xbin, h, 0, 1000);
            plot(x, f{1,2}, 'Color', col.vars.(param)(idx_quant,:), 'LineWidth', 3)

            %[f, x] = ksdensity(bouts_proc.durations_s(quantiles == idx_param), 'Function','pdf', 'NumPoints', 500);
            %length(find(bouts_proc.durations_s(quantiles == idx_param)))

            %plot(x, f, 'Color', col.vars.(param)(1 + idx_param,:), 'LineWidth', 3)
            ylabel('kde');
            ylim([0 1]);

        elseif strcmp(type, 'whisk')
        end
    end
    if strcmp(type, 'whisk')
        boxplot(bouts_proc.durations_s, quantiles)
        ylim([0 10.5]);

    end



% Separate Graph Functions
yticks([0, 0.5, 1])

switch moment
    case 'loom'
        xlim([0 10.5])
        ylim([-0.02 1.02]);
    case 'bsl'
        xlim([-0.05 1.05])
        ylim([0.5 1.0]);
end

if strcmp(type, 'whisk')
    ylim([0 100]);
    xlim([0 num_quantiles + 1]);
    ax.YAxis.Scale = 'log';
end
ax.XAxis.FontSize = 28;
ax.YAxis.FontSize = 28;
set(gca,'linewidth', 3, 'TickDir','out'); % The only other option is 'in'
set(gca,'box','off')
set(gca,'TickLength',[0.02, 0.02])
xlabel('Freeze Duration (s)');
colormap(col.vars.(param)(2:end,:));
clim([0.5 idx_quant + 0.5])

if strcmp(param, 'avg_fs_1s_norm') || strcmp(param, 'avg_sm_freeze_norm') || strcmp(param, 'onsets_loomaligned_norm')

    cell_vec = arrayfun(@(x) ['Q', num2str(x)], 1:idx_quant, 'UniformOutput', false);
    colorbar(ax, 'northoutside', 'Ticks', 1:num_quantiles, 'TickLabels', cell_vec, 'FontSize', 18, 'LineWidth', 2, 'TickLength', 0.1, 'TickDirection', 'none');

elseif strcmp(param, 'nloom')
    cell_vec = {'1-5', '6-10', ; '11-15', '16-20'};
    colorbar(ax, 'northoutside', 'Ticks', 1:num_quantiles, 'TickLabels', cell_vec, 'FontSize', 18, 'LineWidth', 2, 'TickLength', 0.1, 'TickDirection', 'none');

elseif strcmp(param, 'sloom')
    cell_vec = {'Slow', 'Fast'};
    colorbar(ax, 'northoutside', 'Ticks', 1:num_quantiles, 'TickLabels', cell_vec, 'FontSize', 18, 'LineWidth', 2, 'TickLength', 0.1, 'TickDirection', 'none');

elseif strcmp(param, 'moving_flies')
    cell_vec = arrayfun(@(x) [num2str(x)], 0:idx_quant-1, 'UniformOutput', false);
    colorbar(ax, 'northoutside', 'Ticks', 1:num_quantiles, 'TickLabels', cell_vec, 'FontSize', 18, 'LineWidth', 2, 'TickLength', 0.1, 'TickDirection', 'none');
end

text(ax.XLim(2)/2, 0.987, param_str, 'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize', 20, 'Color', col.vars.(param)(idx_quant,:), 'Clipping', 'off');
disp('-------')

if export
    figure_title = sprintf('fd_%s_%s', param, moment);
    exporter(fh, paths, [figure_title, '.png'])
    exporter(fh, paths, [figure_title, '.eps'])
end
end

end