%% Script Name:
% fd_distr_with_param_new
% last update: 18.02

%% Description:
% Plots the distribution of freeze duration as a function of a specific
% experimental parameter

%% Load Data and Preprocess

function fd_distr_withparam_new(varargin)

opt = inputParser;
addParameter(opt, 'bouts', []);
addParameter(opt, 'export', false);
addParameter(opt, 'period', []);
addParameter(opt, 'paths', 'descr/fd_distr');
addParameter(opt, 'params', {'avg_sm_freeze_norm', 'avg_fs_1s_norm', 'sloom', 'nloom'});
addParameter(opt, 'type', 'ecdf');

parse(opt, varargin{:});

bouts_proc = opt.Results.bouts;
export = opt.Results.export;
period = opt.Results.period;
paths = opt.Results.paths;
params = opt.Results.params;
type = opt.Results.type;

if isempty(period)
    if unique(bouts_proc.period) == 1
        period = 'loom';
    elseif unique(bouts_proc.period) == 0
        period = 'bsl';
    end
end

% Load Colors
extra.quantiles = 5;
col = cmapper([], extra.quantiles);

[g, ~] = findgroups(bouts_proc.fly);
c = splitapply(@(x1){cumsum(x1)}, ones(height(bouts_proc),1), g);
bouts_proc.n_freezes = vertcat(c{:});

%%%%%%%%%%%%%%%%%%%%%
% Choose Param here %
%%%%%%%%%%%%%%%%%%%%%

for idx_param = params

    param = idx_param{1};

    if strcmp(param, 'sloom')

    end

    
    if strcmp(param, 'nloom')
        num_quantiles = 4;
        param_str = 'Loom Number';
        thresholds = 0:5:21;

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

    % Load Colors
    extra.quantiles = num_quantiles;
    col = cmapper([], extra.quantiles);

    % Create Figure
    fh = figure('color','w','Position',[100, 100, 400, 500]);
    tl = tiledlayout(10, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile(4, [7 1])
    ax = gca;
    hold on

    for idx_quant = 1:num_quantiles
        if strcmp(type, 'ecdf')
            
            [f, x] = ecdf(bouts_proc.durations_s(quantiles == idx_quant));%, 'censoring', bouts_proc.bout_with_loom(quantiles == idx_quant) > 0);

            plot(x, f, 'Color', col.vars.(param)(1 + idx_quant,:), 'LineWidth', 3)
            ylabel('ecdf');

            switch period
                case 'loom'
                    xlim([0 10.5])
                    ylim([-0.02 1.02]);
                case 'bsl'
                    xlim([-0.05 1.05])
                    ylim([0.5 1.0]);
            end
            text(ax.XLim(2)/2, 0.99, param_str, 'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize', 20, 'Color', col.vars.(param)(idx_quant,:));

        elseif strcmp(type, 'kde')
            
            % RTs{1,1} = bouts_proc.durations_s(quantiles == idx_quant); RTs{2,1} = bouts_proc.durations_s(quantiles == idx_quant);
            % [f] = kreg_single(RTs, RTs, x, xbin, h, 0, 1000);
            % plot(x, f{1,2}, 'Color', col.vars.(param)(idx_quant,:), 'LineWidth', 3)

            %[f, x] = ksdensity(bouts_proc.durations_s(quantiles == idx_param), 'Function','pdf', 'NumPoints', 500);
            %length(find(bouts_proc.durations_s(quantiles == idx_param)))

            %plot(x, f, 'Color', col.vars.(param)(1 + idx_param,:), 'LineWidth', 3)

            histogram(bouts_proc.durations_s(quantiles == idx_quant), min(bouts_proc.durations_s):5/60:max(bouts_proc.durations_s), 'DisplayStyle', 'stairs',...
                'EdgeColor', col.vars.(param)(1 + idx_quant,:), 'LineWidth', 1.2, 'Normalization', 'pdf')

            ylabel('pdf');

            switch period
                case 'loom'
                    xlim([0 10.5])
                    ylim([-0.02 1.52]);
                case 'bsl'
                    xlim([-0.05 1.05])
                    ylim([-0.02 1.52]);
            end

            text(ax.XLim(2)/2, 1.52, param_str, 'HorizontalAlignment','center','VerticalAlignment','top','FontSize', 20, 'Color', col.vars.(param)(idx_quant,:));

        end
    end

    if strcmp(type, 'whisk')
        bh = boxchart(quantiles, bouts_proc.durations_s, 'GroupByColor', quantiles, 'Orientation', 'horizontal', 'BoxWidth', 2.2);
        colororder(col.vars.(param)(2:end, :))
        yticks(1:num_quantiles)
        switch period
            case 'loom'
                xlim([0 10.5])
                ylim([0 5]);
            case 'bsl'
                xlim([-0.05 1.05])
                ylim([0 5]);
        end
    end

    % Separate Graph Functions
   
    ax = gca;
    apply_generic(gca, 'tick_length', 0.02)
    
    xlabel('Freeze Duration (s)');
    colormap(col.vars.(param)(2:end,:));
    clim([0.5 idx_quant + 0.5])

    if strcmp(param, 'avg_fs_1s_norm') || strcmp(param, 'avg_sm_freeze_norm') || strcmp(param, 'onsets_loomaligned_norm')

        cell_vec = arrayfun(@(x) ['Q', num2str(x)], 1:idx_quant, 'UniformOutput', false);
        colorbar(ax, 'northoutside', 'Ticks', 1:num_quantiles, 'TickLabels', cell_vec, 'FontSize', 18, 'LineWidth', 2, 'TickLength', 0.1, 'TickDirection', 'none');

    elseif strcmp(param, 'nloom')
        cell_vec = {'1-5', '6-10', '11-15', '16-20'};
        colorbar(ax, 'northoutside', 'Ticks', 1:num_quantiles, 'TickLabels', cell_vec, 'FontSize', 18, 'LineWidth', 2, 'TickLength', 0.1, 'TickDirection', 'none');

    elseif strcmp(param, 'sloom')
        cell_vec = {'Slow', 'Fast'};
        colorbar(ax, 'northoutside', 'Ticks', 1:num_quantiles, 'TickLabels', cell_vec, 'FontSize', 18, 'LineWidth', 2, 'TickLength', 0.1, 'TickDirection', 'none');

    elseif strcmp(param, 'moving_flies')
        cell_vec = arrayfun(@(x) [num2str(x)], 0:idx_quant-1, 'UniformOutput', false);
        colorbar(ax, 'northoutside', 'Ticks', 1:num_quantiles, 'TickLabels', cell_vec, 'FontSize', 18, 'LineWidth', 2, 'TickLength', 0.1, 'TickDirection', 'none');
    end

    nexttile(1, [3 1])

    if strcmp(param, 'avg_fs_1s_norm') || strcmp(param, 'avg_sm_freeze_norm') || strcmp(param, 'onsets_loomaligned_norm')

        [f, xi] = ksdensity(bouts_proc.(param), 'Bandwidth', 0.01, 'Kernel','normal', 'NumPoints', 500);
        for idx_quant = 1:num_quantiles
            fill_between(xi(xi >= thresholds(idx_quant) & xi < thresholds(idx_quant + 1)) , zeros(length(xi(xi >= thresholds(idx_quant) & xi < thresholds(idx_quant + 1))), 1), f(xi >= thresholds(idx_quant) & xi < thresholds(idx_quant + 1)),[], 'FaceColor',col.vars.(param)(1 + idx_quant,:),'LineStyle','none')
            % plot(xi,f)
        end
        xlim([0 2])
        xticks(0:3)

    elseif strcmp(param, 'nloom') || strcmp(param, 'sloom') || strcmp(param, 'moving_flies')

        bh = bar(1:num_quantiles, histcounts(quantiles), 'FaceColor', 'flat', 'EdgeColor', 'none');
        bh.CData = col.vars.(param)(2:end,:);
        xticks([])
        xlim([0.5, num_quantiles + 0.5])

    end

    ax = gca;
    set(gca,'box','off')
    ax.YAxis.Visible = 'off';
    ax.XAxis.Visible = 'on';
    ax.XAxis.FontSize = 18;
    set(gca,'linewidth', 2, 'TickDir','out'); % The only other option is 'in'
    set(gca,'box','off')

    if export
    figure_title = sprintf('fd_%s_%s', param, moment);
    exporter(fh, paths, [figure_title, '.pdf'])
    end
    
end

end

