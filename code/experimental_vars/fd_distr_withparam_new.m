function fh = fd_distr_withparam_new(varargin)

%% ===================== Parse inputs =====================
opt = inputParser;
addParameter(opt,'bouts',[]);
addParameter(opt,'export',false);
addParameter(opt,'period',[]);
addParameter(opt,'paths','descr/fd_distr');
addParameter(opt,'params',{'avg_sm_freeze_norm','avg_fs_1s_norm','sloom','nloom'});
addParameter(opt,'type','ecdf');
addParameter(opt,'check_quantiles',false);
addParameter(opt,'sloom_to_plot','all');

parse(opt,varargin{:});

bouts_proc      = opt.Results.bouts;
export          = opt.Results.export;
period          = opt.Results.period;
paths           = opt.Results.paths;
params          = opt.Results.params;
type            = opt.Results.type;
check_quantiles = opt.Results.check_quantiles;
sloom_to_plot   = opt.Results.sloom_to_plot;

%% ===================== Infer period =====================
if isempty(period)
    u = unique(bouts_proc.period);
    if numel(u) ~= 1
        error('Mixture of baseline and loom periods');
    end
    if u == 1
        period = 'loom';
    else
        period = 'bsl';
    end
end

%% ===================== Determine which sloom values to plot =====================
all_sloom = unique(bouts_proc.ls)';

switch sloom_to_plot
    case 'all'
        sloom_vals = all_sloom;
    case 'slow'
        sloom_vals = min(all_sloom);
    case 'fast'
        sloom_vals = max(all_sloom);
    otherwise
        sloom_vals = sloom_to_plot;
end

n_sloom = numel(sloom_vals);

%% ===================== Loop over parameters =====================
for idx_param = params
    param = idx_param{1};
    is_sloom_param = strcmp(param,'sloom');

    [num_quantiles, thresholds] = param_settings(param, bouts_proc);
    quantiles = discretize(bouts_proc.(param), thresholds);

    %% ===================== Colormap =====================
    col  = cmapper([], num_quantiles);
    cmap = get_var_cmap(col, param);

    %% ===================== Figure =====================
    if n_sloom == 1 || is_sloom_param
        fig_width = 380;
    else
        fig_width = 500;
    end

    fh = figure('Color','w','Position',[100 100 fig_width 500]);
    tiledlayout(10,2,'TileSpacing','tight','Padding','compact');

    %% ===================== Control distribution =====================
    nexttile(1,[3 2]); hold on
    plot_control_distribution( ...
        bouts_proc.(param), quantiles, thresholds, cmap, param)

    ax = gca;
    pad_ylim(ax,0.075)
    apply_generic(ax,'no_y',true,'no_x',false,'font_size',20)

    mid_q = ceil(num_quantiles/2) + 1;

    if strcmp(param,'avg_sm_freeze_norm')
        xlabel('Avg. Social Motion','Interpreter','none','Color',cmap(mid_q,:))
    elseif strcmp(param,'avg_fs_1s_norm')
        xlabel('Focal Speed before Loom','Interpreter','none','Color',cmap(mid_q,:))
    elseif strcmp(param, 'sloom')
        xlabel('Loom Speed','Interpreter','none','Color',cmap(mid_q,:))
    end

    colormap(cmap(2:end,:));
    clim([0 num_quantiles])
    add_quantile_colorbar(ax,param,num_quantiles)

    ax.XAxis.Exponent = 0;

    %% ===================== Duration distributions =====================
    ax_bottom = [];

    if is_sloom_param
        % ===== slow & fast loom in SAME tile =====
        nexttile([7 2]); hold on

        for idx_sloom = sloom_vals
            bouts_sloom = bouts_proc(bouts_proc.ls==idx_sloom,:);
            qmask = quantiles(bouts_proc.ls==idx_sloom);

            plot_duration_distribution( ...
                bouts_sloom.durations_s, qmask, ...
                cmap, num_quantiles, type, period);
        end

        xlabel('Freeze Duration (s)');
        ylabel(type);
        apply_generic(gca,'tick_length',0.025,'font_size', 26, ...
            'yticks',[0 1],'xticks',[0 10])

        ax_bottom = gca;

    else
        % ===== one tile per loom speed =====
        if n_sloom == 1
            tile_width = 2;
        else
            tile_width = 1;
        end

        for idx_sloom = sloom_vals
            bouts_sloom = bouts_proc(bouts_proc.ls==idx_sloom,:);
            qmask = quantiles(bouts_proc.ls==idx_sloom);
            censored = bouts_sloom.censored_contacts;
            nexttile([7 tile_width]); hold on
            plot_duration_distribution( ...
                bouts_sloom.durations_s, qmask, ...
                cmap, num_quantiles, type, period, censored);

            ylabel(type);
            apply_generic(gca,'tick_length',0.025,'font_size',26, ...
                'yticks',[0 1],'xticks',[0 10])

            if idx_sloom == min(all_sloom) && strcmp(sloom_to_plot, 'all')
                txt = 'Slow Loom';
            else
                txt = 'Fast Loom';
            end

            text(gca,1.015,0.015,txt, ...
                'Units','normalized', ...
                'FontSize',16, ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','bottom');

            ax_bottom(end+1) = gca; %#ok<AGROW>
        end
        xlabel('Freeze Duration (s)');

    end

    if numel(ax_bottom) > 1
        linkaxes(ax_bottom,'y')
    end

    %% ===================== Export =====================
    if export
        exporter(fh, paths, sprintf('fd_%s_%s.pdf', param, period))
                exporter(fh, paths, sprintf('fd_%s_%s.png', param, period))

    end
end
end

%% =====================================================================
%% ============================ HELPERS ================================
%% =====================================================================

function [nq, thresholds] = param_settings(param,T)
switch param
    case 'nloom'
        nq = 4; thresholds = 0:5:21;
    case 'sloom'
        nq = 2; thresholds = 12.5:25:75;
    case 'moving_flies'
        nq = 5; thresholds = -0.5:1:5.5;
    otherwise
        nq = 4;
        thresholds = quantile(T.(param),linspace(0,1,nq+1));
end
end

% ---------------------------------------------------------------------

function plot_control_distribution(data,quantiles,thresholds,cmap,param)

if ismember(param,{'nloom','sloom','moving_flies'})
    bh = bar(1:max(quantiles),histcounts(quantiles), ...
        'FaceColor','flat','EdgeColor','none');
    bh.CData = cmap(2:end,:);
    xlim([.5 max(quantiles)+.5])
    xticks([])
else
    iqr_val = iqr(data);
    if iqr_val > 0
        bw = 2 * iqr_val / numel(data)^(1/2.1);
        nbins = max(20, ceil(range(data)/bw));
    else
        nbins = 50;
    end

    [c,e] = histcounts(data,nbins);
    centers = e(1:end-1)+diff(e)/2;

    b = bar(centers,c,1,'FaceColor','flat','EdgeColor','none');
    for q = 1:max(quantiles)
        mask = centers>=thresholds(q) & centers<thresholds(q+1);
        b.CData(mask,:) = repmat(cmap(1+q,:),sum(mask),1);
    end

    lo = prctile(data,1);
    hi = prctile(data,99);
    if lo < hi
        xlim([lo hi])
    end
end
end

% ---------------------------------------------------------------------

function plot_duration_distribution(durations,qmask,cmap,nq,type,period, censored)

switch type
    case 'cumulative'
        for q = 1:nq
            d = durations(qmask==q);
            c = censored(qmask==q);
            if isempty(d), continue, end
            sum(c)
            [f,x] = ecdf(d, 'Censoring', c);
            plot(x, f,'LineWidth', 4 ,'Color',cmap(1+q,:) )
        end
        ylim([0 1])
        pad_ylim(gca,0.025)

    case 'kde'
        edges = min(durations):5/60:max(durations);
        for q = 1:nq
            d = durations(qmask==q);
            if isempty(d), continue, end
            histogram(d,edges,'Normalization','pdf', ...
                'DisplayStyle','stairs', ...
                'LineWidth',1.5, ...
                'EdgeColor',cmap(1+q,:));
        end
        pad_ylim(gca,0.005)

    case 'whisk'
        for q = 1:nq
            d = durations(qmask==q);
            if isempty(d), continue, end
            boxchart(q*ones(size(d)),d, ...
                'Orientation','horizontal', ...
                'BoxWidth',0.5, ...
                'BoxFaceColor',cmap(1+q,:), ...
                'MarkerStyle','none', ...
                'LineWidth',1.2);
        end
        yticks(1:nq)
        yticklabels(arrayfun(@(x)['Q' num2str(x)],1:nq,'uni',0))
        ylim([0.5 nq+0.5])
end

switch period
    case 'loom', xlim([0 10.5])
    case 'bsl',  xlim([-0.05 1.05])
end
end

% ---------------------------------------------------------------------

function cmap = get_var_cmap(col,varname)
if isfield(col.vars,varname)
    cmap = col.vars.(varname);
else
    cmap = col.vars.default;
end
end

% ---------------------------------------------------------------------

function add_quantile_colorbar(ax,param,nq)
switch param
    case 'nloom'
        labels = {'1–5','6–10','11–15','16–20'};
    case 'sloom'
        labels = {'Slow','Fast'};
    otherwise
        labels = arrayfun(@(x)['Q' num2str(x)],1:nq,'uni',0);
end

colorbar(ax,'southoutside', ...
    'Ticks',0.5:nq+0.5, ...
    'TickLabels',labels, ...
    'Limits',[0 nq], ...
    'TickDirection','none', ...
    'FontSize', 22, ...
    'LineWidth',2)
end

% ---------------------------------------------------------------------

function pad_ylim(ax,fraction)
yl = ax.YLim;
p = fraction * diff(yl);
ax.YLim = [yl(1)-p yl(2)+p];
end
