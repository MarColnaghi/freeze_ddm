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

parse(opt,varargin{:});

bouts_proc = opt.Results.bouts;
export = opt.Results.export;
period = opt.Results.period;
paths = opt.Results.paths;
params = opt.Results.params;
type = opt.Results.type;
check_quantiles = opt.Results.check_quantiles;

%% ===================== Infer period =====================
if isempty(period)
    u = unique(bouts_proc.period);
    if numel(u) ~= 1
        error('Mixture of baseline and loom periods');
    end
    period = ternary(u==1,'loom','bsl');
end

%% ===================== Loop over parameters =====================
for idx_param = params
    param = idx_param{1};

    [num_quantiles, thresholds] = param_settings(param, bouts_proc);
    quantiles = discretize(bouts_proc.(param), thresholds);

    if check_quantiles && ismember('sloom_norm', bouts_proc.Properties.VariableNames)

        fprintf('\nParam: %s\n', param);

        u_sloom = unique(bouts_proc.sloom_norm);

        thr_by_sloom = nan(numel(u_sloom), num_quantiles+1);

        for i = 1:numel(u_sloom)
            s = u_sloom(i);
            idx = bouts_proc.sloom_norm == s;

            data_s = bouts_proc.(param)(idx);

            if numel(data_s) < num_quantiles
                warning('Not enough data for sloom = %.2f', s);
                continue
            end

            thr_by_sloom(i,:) = quantile( ...
                data_s, linspace(0,1,num_quantiles+1));

            fprintf('  sloom = %.2f | ', s);
            fprintf('%.3f ', thr_by_sloom(i,:));
            fprintf('\n');
        end

        % ---- Optional plot (very helpful)
        figure('Color','w'); hold on
        for i = 1:numel(u_sloom)
            plot(0:num_quantiles, thr_by_sloom(i,:), '-o', ...
                'DisplayName', sprintf('sloom = %.2f', u_sloom(i)));
        end
        xlabel('Quantile index')
        ylabel(strrep(param,'_',' '))
        title(sprintf('Quantile thresholds by loom speed: %s', param), ...
            'Interpreter','none')
        legend('Location','best')
    end

    col  = cmapper([], num_quantiles);
    cmap = get_var_cmap(col, param);   % SAFE ACCESS

    %% ===================== Figure =====================
    fh = figure('Color','w','Position',[100 100 700 500]);
    tiledlayout(10,2,'TileSpacing','tight','Padding','compact');

    %% ===================== Tile 1: control variable =====================
    nexttile(1,[3 2]); hold on
    plot_control_distribution( ...
        bouts_proc.(param), quantiles, thresholds, cmap, param)

    ax = gca;
    pad_ylim(ax,0.075)
    apply_generic(ax,'no_y',true,'no_x',false,'font_size',18)

    colormap(cmap(2:end,:));
    clim([0 num_quantiles])
    add_quantile_colorbar(ax,param,num_quantiles)

    ax.XAxis.Exponent = 0;   % disable ×10^n
    ax.XTickLabelRotation = 0;

    %% ===================== Duration distributions =====================
    ax_bottom = gobjects(0);

    for idx_sloom = unique(bouts_proc.sloom_norm)'

        bouts_sloom = bouts_proc(bouts_proc.sloom_norm==idx_sloom,:);
        qmask = quantiles(bouts_proc.sloom_norm==idx_sloom);

        nexttile([7 1]); hold on
        plot_duration_distribution( ...
            bouts_sloom.durations_s, qmask, ...
            cmap, num_quantiles, type, period);

        xlabel('Freeze Duration (s)');
        ylabel(type);
        apply_generic(gca,'tick_length',0.02,'font_size',20)

        % ---- Condition label
        txt = ternary(idx_sloom==1,'Slow Loom','Fast Loom');
        text(gca,1.015,0.015,txt, ...
            'Units','normalized', ...
            'FontSize',16, ...
            'HorizontalAlignment','right', ...
            'VerticalAlignment','bottom');

        ax_bottom(end+1) = gca;
    end

    if numel(ax_bottom)>1
        linkaxes(ax_bottom,'y')
    end

    mid_q = ceil(num_quantiles/2) + 1;  % middle quantile index in the colormap
    text(ax, 0.98, 0.98, strrep(param,'_',' '), ...
        'Units','normalized', ...
        'FontSize', 16, ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','top', ...
        'Color', cmap(mid_q,:), ...  % color of middle quantile
        'Interpreter','none');   % prevents underscores from being subscripts

    %% ===================== Export =====================
    if export
        exporter(fh, paths, sprintf('fd_%s_%s.pdf', param, period))
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
    xlim([.5 max(quantiles) + .5])
    xticks([])
else
    % ---- Adaptive binning (Freedman–Diaconis)
    iqr_val = iqr(data);
    if iqr_val > 0
        bw = 2 * iqr_val / numel(data)^(1/3);
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


% ---- Robust x-limits for long tails
lo = prctile(data,1);
hi = prctile(data,99);
if lo < hi
    xlim([lo hi])
end
end
end
% ---------------------------------------------------------------------

function plot_duration_distribution(durations,qmask,cmap,nq,type,period)

switch type
    case 'ecdf'
        for q = 1:nq
            d = durations(qmask==q);
            if isempty(d), continue, end
            [f,x] = ecdf(d);
            plot(x,f,'LineWidth',3,'Color',cmap(1+q,:))
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
    'FontSize',18, ...
    'LineWidth',2)
end

% ---------------------------------------------------------------------

function pad_ylim(ax,fraction)
yl = ax.YLim;
p = fraction*diff(yl);
ax.YLim = [yl(1)-p yl(2)+p];
end

% ---------------------------------------------------------------------

function out = ternary(cond,a,b)
if cond, out = a; else, out = b; end
end
