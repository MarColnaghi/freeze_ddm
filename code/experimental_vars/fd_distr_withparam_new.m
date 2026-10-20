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
addParameter(opt,'scale','linear');
addParameter(opt,'control_var',true);

parse(opt,varargin{:});

bouts_proc      = opt.Results.bouts;
export          = opt.Results.export;
period          = opt.Results.period;
paths           = opt.Results.paths;
params          = opt.Results.params;
type            = opt.Results.type;
check_quantiles = opt.Results.check_quantiles;
sloom_to_plot   = opt.Results.sloom_to_plot;
scale           = opt.Results.scale;
control_var     = opt.Results.control_var;
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
    if control_var
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

    end
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
            if strcmp(type, 'cumulative')

                xl = xlim(gca); yl = ylim(gca);
                text(gca, xl(2) - 0.02*diff(xl), yl(1) + 0.03*diff(yl), txt, ...
                    'FontSize',16, ...
                    'HorizontalAlignment','right', ...
                    'VerticalAlignment','bottom');

                ax_bottom(end+1) = gca; %#ok<AGROW>

            elseif strcmp(type, 'kde')
                xl = xlim(gca); yl = ylim(gca);
                text(gca, xl(2) - 0.02*diff(xl), yl(2) - 0.1*diff(yl), txt, ...
                    'FontSize',16, ...
                    'HorizontalAlignment','right', ...
                    'VerticalAlignment','bottom');

                ax_bottom(end+1) = gca; %#ok<AGROW>
            end
        end
        xlabel(tl,'Freeze Duration (s)','FontSize',26);

    end

    if numel(ax_bottom) > 1
        linkaxes(ax_bottom,'y')
    end

    % scale 'log': the cumulative has a log x axis and the kde a log y axis. The
    % cumulative gets reference lines at 0.2, 0.5 and 1 s, labelled when the axis is
    % log, bare when linear. The kde never has them
    log_x = strcmp(scale, 'log') && strcmp(type, 'cumulative');
    ax_duration = findobj(fh, 'Type', 'axes')';
    if control_var
        ax_duration = ax_duration(1:end-1); % findobj lists the newest first, the control tile is the first one made
    end
    for idx_ax = 1:numel(ax_duration)
        ax = ax_duration(idx_ax);
        if strcmp(scale, 'log') && ~log_x
            ax.YScale = 'log';
        end
        if log_x
            ax.XScale = 'log';
            xlim(ax, [0.1 10])
        else
            xlim(ax, [0 10])
        end

        if strcmp(type, 'cumulative')
            if strcmp(scale, 'log')
                xline(ax, [0.2, 0.5, 1], 'k-', {'0.2s', '0.5s', '1s'}, ...
                    'LabelHorizontalAlignment', 'right', ...
                    'LabelVerticalAlignment', 'top', 'FontSize', 16, 'LineWidth',2);
            else
                xline(ax, [0.2, 0.5, 1], 'k-', 'LineWidth',2);
            end
        end
    end

    %% ===================== Export =====================
    if export 
        if strcmp(scale, 'log')
        exporter(fh, paths, sprintf('fd_%s_%s_log.pdf', param, period))
                exporter(fh, paths, sprintf('fd_%s_%s.png', param, period))
        else
            exporter(fh, paths, sprintf('fd_%s_%s.pdf', param, period))
                exporter(fh, paths, sprintf('fd_%s_%s.png', param, period))
        end
    end
end
end


%% =====================================================================
%% ============================ HELPERS ================================
%% =====================================================================

function [nq, thresholds] = param_settings(param,T)
switch param
    case {'nloom', 'ln'}
        % Edges between the loom numbers, so that 1-5, 6-10, 11-15 and 16-20 are the bins.
        % ln is nloom / 10
        nq = 4; thresholds = 0.5:5:20.5;
        if strcmp(param, 'ln'), thresholds = thresholds / 10; end
    case {'ls','jump_start','late_onset'}
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
    case 'late_onset',         s = 'Onset after split';
    case {'nloom', 'ln'},      s = 'Loom Number';
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

% Discrete variables: one bar per level, the count of freezes, under the colorbar
% labels. The others get the histogram of the variable with the quantiles coloured
if ismember(param,{'ln','ls','nloom','moving_flies','jump_start','late_onset'})
    bh = bar(1:max(quantiles),histcounts(quantiles,0.5:1:max(quantiles)+0.5), ...
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
            stairs(x, f,'LineWidth', 4 ,'Color',cmap(1+q,:) )
        end
        ylim([0 1])
        pad_ylim(gca,0.025)

    case 'kde'
        % The Kaplan-Meier mass of every bin, so a freeze that was cut counts where it
        % was seen to last and not as if it ended there. A bin that reaches past the
        % longest observation, when that one was cut, is left blank, and a curve adds
        % up to the F(t) of the cumulative plot, not to 1
                bin_w = 10 / 60;
        n_bin = max(1, ceil((max(durations) - min(durations)) / bin_w));
        edges = min(durations) + (0:n_bin) * bin_w;
        for q = 1:nq
            d = durations(qmask == q);
            c = censored(qmask == q);
            if isempty(d) || all(c), continue, end
            [mass, covered] = km_bin_mass(d, c, edges);
            h = mass / bin_w;
            h(~covered) = NaN;
            stairs(edges, [h h(end)], 'LineWidth', 1.25, 'Color', cmap(1 + q, :))
        end
        pad_ylim(gca, 0.005)
        ylim([-0.05 2.05])

    case 'whisk'
        % Does not take the censoring into account
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

function [mass,covered] = km_bin_mass(d,c,edges)
% Kaplan-Meier mass of every bin: the steps of the curve that fall in it. covered is
% false for the bins that reach past the longest observation when that one was cut,
% their mass is not known
covered = true(1,numel(edges)-1);
[f,x] = ecdf(d,'Censoring',c);
step = diff(f);
b = discretize(x(2:end),edges);
ok = ~isnan(b);
mass = accumarray(b(ok),step(ok),[numel(edges)-1,1])';
if f(end) < 1-1e-9
    covered = edges(2:end) <= max(d)+1e-9;
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
    case {'nloom', 'ln'}
        % The first bin starts at the first loom in the data, 2 when loom 1 is left out
        labels = {sprintf('%d–5',max(1,min(T.nloom))),'6–10','11–15','16–20'};
    case {'sloom','ls'}
        labels = {'Slow','Fast'};
    case 'late_onset'
        labels = {'Early','Late'};
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
