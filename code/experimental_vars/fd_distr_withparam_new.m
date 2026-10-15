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
    is_sloom_param = strcmp(param, 'ls');

    [num_quantiles, thresholds] = param_settings(param, bouts_proc);
    quantiles = discretize(bouts_proc.(param), thresholds);

    %% ===================== Colormap =====================
    col  = cmapper([], num_quantiles);
    cmap = get_var_cmap(col, param);

    %% ===================== Figure =====================
    if n_sloom == 1 || is_sloom_param
        fig_width = 380;
    else
        fig_width = 600;
    end

    fh = figure('Color','w','Position',[100 100 fig_width 500]);
    tl = tiledlayout(10,2,'TileSpacing','tight','Padding','compact');

    %% ===================== Control distribution =====================
    nexttile(1,[3 2]); hold on
   plot_control_distribution( ...
       bouts_proc.(param), quantiles, thresholds, cmap, param)

    ax = gca;
    pad_ylim(ax,0.075)
    apply_generic(ax,'no_y',true,'no_x',false,'font_size',20)

    mid_q = ceil(num_quantiles/2) + 1;

    xlabel(control_label(param),'Interpreter','none','Color',cmap(mid_q,:))

    colormap(cmap(2:end,:));
    clim([0 num_quantiles])
    add_quantile_colorbar(ax,param,num_quantiles,bouts_proc)

    ax.XAxis.Exponent = 0;

    %% ===================== Duration distributions =====================
    ax_bottom = [];

    if is_sloom_param
        % ===== slow & fast loom in SAME tile =====
        nexttile([7 2]); hold on

        for idx_sloom = sloom_vals
            bouts_sloom = bouts_proc(bouts_proc.ls==idx_sloom,:);
            qmask = quantiles(bouts_proc.ls==idx_sloom);
            censored = bouts_sloom.is_censored;

            plot_duration_distribution( ...
                bouts_sloom.durations_s, qmask, ...
                cmap, num_quantiles, type, period, censored);
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
            censored = bouts_sloom.is_censored;
            nexttile([7 tile_width]); hold on
            plot_duration_distribution( ...
                bouts_sloom.durations_s, qmask, ...
                cmap, num_quantiles, type, period, censored);

            ylabel(type);
            apply_generic(gca,'tick_length',0.025,'font_size',26, ...
                'yticks',[0 1],'xticks',[0 10])

            % ls is 0 for the slow loom and 1 for the fast one. Going by the smallest
            % value in the data called the only speed of a 50 cm dataset slow
            if idx_sloom == 0
                txt = 'Slow Loom';
            else
                txt = 'Fast Loom';
            end

            % Data coordinates, which render the same everywhere, unlike normalized
            % text in a tiled layout
            xl = xlim(gca); yl = ylim(gca);
            text(gca, xl(2) - 0.02*diff(xl), yl(1) + 0.03*diff(yl), txt, ...
                'FontSize',16, ...
                'HorizontalAlignment','right', ...
                'VerticalAlignment','bottom');

            ax_bottom(end+1) = gca; %#ok<AGROW>
        end
        xlabel(tl,'Freeze Duration (s)','FontSize',26);

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
        % Edges between the loom numbers, so that 1-5, 6-10, 11-15 and 16-20 are the bins
        nq = 4; thresholds = 0.5:5:20.5;
    case {'ls','jump_start'}
        % binary: quantile() edges collapse to [0 0 0 0 1] and discretize
        % puts every bout in the top bin
        nq = 2; thresholds = -0.5:1:1.5;
    case 'moving_flies'
        nq = 5; thresholds = -0.5:1:5.5;
    otherwise
        nq = 4;
        thresholds = quantile(T.(param),linspace(0,1,nq+1));
end
end

% ---------------------------------------------------------------------

function s = control_label(param)
% What the control variable is called on the figure; the column name if it has no entry
switch param
    case 'avg_sm_freeze_norm', s = 'Avg. Social Motion';
    case 'avg_fs_1s_norm',     s = 'Focal Speed before Loom';
    case 'ls',                 s = 'Loom Speed';
    case 'moving_flies',       s = 'Moving Flies';
    case 'nloom',              s = 'Loom Number';
    case 'sm_1s',              s = 'Social Motion (1st second)';
    case 'sm_2s',              s = 'Social Motion (1st 2 s)';
    case 'sm',                 s = 'Social Motion (freeze)';
    case 'avg_fs_1s',          s = 'Focal Speed (1 s before freeze)';
    case 'avg_fs_loom',        s = 'Focal Speed (before previous loom)';
    otherwise,                 s = strrep(param, '_', ' ');
end
end

% ---------------------------------------------------------------------

function plot_control_distribution(data,quantiles,thresholds,cmap,param)

if ismember(param,{'ln','ls','moving_flies', 'jump_start'})
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
            [f,x] = ecdf(d, 'Censoring', c);
            plot(x, f,'LineWidth', 4 ,'Color',cmap(1+q,:) )
        end
        ylim([0 1])
        pad_ylim(gca,0.025)

    case 'kde'
        edges = min(durations):20/60:max(durations);
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
    case 'loom', xlim([0 10])
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

function add_quantile_colorbar(ax,param,nq,T)
switch param
    case 'nloom'
        % The first bin starts at the first loom in the data, 2 when loom 1 is left out
        labels = {sprintf('%d–5',max(1,min(T.nloom))),'6–10','11–15','16–20'};
    case {'sloom','ls'}
        labels = {'Slow','Fast'};
    case 'moving_flies'
        labels = arrayfun(@num2str,0:nq-1,'uni',0);
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
