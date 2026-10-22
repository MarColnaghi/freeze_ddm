function fh = freeze_break_sm_by_param(varargin)
% Social motion the focal fly experiences, aligned to the moment it breaks its
% freeze (t = 0, the first frame after the freeze), split by the levels of a
% control variable: one mean curve with a SEM band for each level, as
% fd_distr_withparam_new does for the freeze duration.
%
%   fh = freeze_break_sm_by_param('bouts', bouts_proc, 'params', {'all', 'ls', 'moving_flies'})
%
% bouts is the table of loom_evoked_freezes, whose contact censoring
% (impose_contact_threshold) is already applied. Only the freezes that ended by
% themselves (not is_censored) have a break to align to, so those are the ones
% averaged. Levels are those of fd_distr_withparam_new: the quantile edges are taken on
% the whole table, censored freezes included, so Q1-Q4 are the same freezes as in
% the duration figures.
%
% The signal is read once, after the censoring and only for the freezes that are
% kept, which keeps the matrix at the width of the longest kept freeze.
%
% Options
%   params             {'all'}  columns of bouts to split by, one figure each. 'all'
%                               is a single curve with no split
%   sloom_to_plot      'all'    'all', 'slow', 'fast' or the ls value, one tile per
%                               loom speed unless the parameter is ls itself
%   split_sloom        true     false pools the loom speeds in one tile (the speeds
%                               left by sloom_to_plot)
%   error              'fly'    'fly' averages the freezes of a fly first and takes the
%                               SEM over flies, 'freeze' treats every freeze as
%                               independent (the band is then too narrow: a fly
%                               contributes several freezes)
%   delay_after        90       frames of signal kept after the break
%   xlims_s            [-6 1.5] seconds around the break
%   min_n              20       a lag is drawn when at least this many flies (freezes
%                               for error 'freeze') contribute to it
%   norm_factor        10       social motion is divided by it, as in the other figures
%   mask_post_contact  true     the signal after the break is dropped from the first
%                               frame the nearest fly is closer than
%                               contact_threshold. The contact censoring only looks
%                               at the freeze, and a contact spikes the signal
%   contact_threshold  []       taken from the table (impose_contact_threshold
%                               writes it) and 70 if it is not there
%   export, paths               as in fd_distr_withparam_new

opt = inputParser;
addParameter(opt, 'bouts', []);
addParameter(opt, 'params', {'all'});
addParameter(opt, 'export', false);
addParameter(opt, 'paths', 'social_motion/event_locked');
addParameter(opt, 'sloom_to_plot', 'all');
addParameter(opt, 'split_sloom', true);
addParameter(opt, 'error', 'fly', @(x) ismember(x, {'fly', 'freeze'}));
addParameter(opt, 'delay_after', 90);
addParameter(opt, 'xlims_s', [-6 1.5]);
addParameter(opt, 'min_n', 20);
addParameter(opt, 'norm_factor', 10);
addParameter(opt, 'mask_post_contact', true);
addParameter(opt, 'contact_threshold', []);
parse(opt, varargin{:});

all_bouts         = opt.Results.bouts;
params            = opt.Results.params;
export            = opt.Results.export;
paths             = opt.Results.paths;
sloom_to_plot     = opt.Results.sloom_to_plot;
split_sloom       = opt.Results.split_sloom;
by_fly           = strcmp(opt.Results.error, 'fly');
delay_after       = opt.Results.delay_after;
xlims_s           = opt.Results.xlims_s;
min_n             = opt.Results.min_n;
norm_factor       = opt.Results.norm_factor;
mask_post_contact = opt.Results.mask_post_contact;
contact_threshold = opt.Results.contact_threshold;

camera_fps = 60;

if ~ismember('is_censored', all_bouts.Properties.VariableNames)
    error('bouts has no is_censored: pass the table of loom_evoked_freezes');
end
if isempty(contact_threshold)
    contact_threshold = 70;
    ud = all_bouts.Properties.UserData;
    if isstruct(ud) && isfield(ud, 'contact_threshold')
        contact_threshold = ud.contact_threshold;
    end
end

%% ===================== Break-aligned signal =====================
% Right-aligned: each row is a freeze, the last columns sit at its break plus
% delay_after, NaNs pad the front of the shorter ones
freezes = all_bouts(~all_bouts.is_censored, :);

sm = extract_sm_from_bouts(freezes, 'type', 'onlyfreeze', 'output_type', 'mat', ...
    'align', 'offset', 'delay_after', delay_after, 'norm_factor', norm_factor, ...
    'cache', 'motion_cache');

% The last column is frame ends - 1 + delay_after, and ends is the first frame after
% the freeze, which is t = 0
t_frames = (-size(sm, 2) + 1 : 0) + delay_after - 1;
t_s = t_frames / camera_fps;

if mask_post_contact
    dist = extract_sm_from_bouts(freezes, 'type', 'onlyfreeze', 'output_type', 'mat', ...
        'align', 'offset', 'delay_after', delay_after, 'cache', 'mindist_cache');
    contact_after = (dist < contact_threshold) & (t_frames >= 0);
    sm(cumsum(contact_after, 2) > 0) = NaN;
    fprintf('%d of %d freezes touched within %.1f s of the break, their signal is dropped from then on\n', ...
        sum(any(contact_after, 2)), height(freezes), delay_after / camera_fps);
end

%% ===================== Loom speeds =====================
all_sloom = unique(all_bouts.ls)';
switch sloom_to_plot
    case 'all',  sloom_vals = all_sloom;
    case 'slow', sloom_vals = min(all_sloom);
    case 'fast', sloom_vals = max(all_sloom);
    otherwise,   sloom_vals = sloom_to_plot;
end

%% ===================== One figure per parameter =====================
for idx_param = params
    param = idx_param{1};
    is_ls_param = strcmp(param, 'ls');

    if strcmp(param, 'all')
        num_levels = 1;
        level = ones(height(freezes), 1);
        col = cmapper();
        cmap = repmat(col.var.sm, 2, 1);
        labels = {'All freezes'};
    else
        [num_levels, edges] = param_settings(param, all_bouts);
        level = discretize(freezes.(param), edges);
        col = cmapper([], num_levels);
        cmap = get_var_cmap(col, param);
        labels = level_labels(param, num_levels, all_bouts);
    end

    % ls itself is the split, so its two levels share one tile, and 'all' pools the speeds
    if is_ls_param || strcmp(param, 'all') || ~split_sloom
        panels = {sloom_vals};
    else
        panels = num2cell(sloom_vals);
    end
    n_panels = numel(panels);


    fh = figure('Color', 'w', 'Position', [100 100 150 + 420 * n_panels 350]);
    tl = tiledlayout(1, n_panels, 'TileSpacing', 'compact', 'Padding', 'compact');
    axes_drawn = gobjects(1, n_panels);

    for idx_panel = 1:n_panels
        ax = nexttile; hold(ax, 'on');
        axes_drawn(idx_panel) = ax;
        in_panel = ismember(freezes.ls, panels{idx_panel});

        h = gobjects(1, num_levels);
        entries = cell(1, num_levels);
        for q = 1:num_levels
            rows = in_panel & level == q;
            if ~any(rows), continue, end

            M = sm(rows, :);
            if by_fly
                [~, ~, g] = unique(freezes.fly(rows));
                M = splitapply(@(x) median(x, 1, 'omitnan'), M, g);
            end
            n = sum(~isnan(M), 1);
            avg = median(M, 1, 'omitnan');
            sem = std(M, 0, 1, 'omitnan') ./ sqrt(n);

            % Only the lags that enough flies reach, which cuts the noisy deep past
            valid = n >= min_n & t_s >= xlims_s(1) & t_s <= xlims_s(2);
            if ~any(valid), continue, end

            c = cmap(1 + q, :);
            fill(ax, [t_s(valid), fliplr(t_s(valid))], [avg(valid) + sem(valid), fliplr(avg(valid) - sem(valid))], c, ...
                'EdgeColor', 'none', 'FaceAlpha', 0.2, 'HandleVisibility', 'off');
            h(q) = plot(ax, t_s(valid), avg(valid), 'LineWidth', 3, 'Color', c);
            entries{q} = sprintf('%s (n = %d)', labels{q}, sum(rows));
        end

        xline(ax, 0, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');

        ticks = ceil(xlims_s(1)) : 2 : floor(xlims_s(2));

        tick_labels = arrayfun(@num2str, ticks, 'UniformOutput', false);
        tick_labels(ticks == 0) = {'Break'};
        xticklabels(ax, tick_labels)

        drawn = isgraphics(h);
        if any(drawn)
            legend(ax, h(drawn), entries(drawn), 'Location', 'northwest', 'Box', 'off', 'FontSize', 14);
        end

        if n_panels > 1
            if panels{idx_panel} == 0, txt = 'Slow Loom'; else, txt = 'Fast Loom'; end
            title(ax, txt, 'FontWeight', 'normal', 'FontSize', 18)
        end
        if idx_panel == 1
            ylabel(ax, 'Social motion')
        end
    end

    if n_panels > 1
        linkaxes(axes_drawn, 'y')
    end
    xlabel(ax, 'Time (s)')
    ylim([0 .6])
    apply_generic(ax, 'font_size', 24, 'tick_length', 0.02, 'line_width', 2.5, ...
        'xlims', xlims_s, 'xticks', ticks);

    if export
        exporter(fh, paths, sprintf('freeze_break_sm_%s.pdf', param))
        exporter(fh, paths, sprintf('freeze_break_sm_%s.png', param))
    end
end
end


%% =====================================================================
%% ============================ HELPERS ================================
%% =====================================================================
% The levels of a parameter, as in fd_distr_withparam_new, so that the same
% quantiles and colours come out

function [nq, thresholds] = param_settings(param, T)
switch param
    case {'nloom', 'ln'}
        nq = 4; thresholds = 0.5:5:20.5;
        if strcmp(param, 'ln'), thresholds = thresholds / 10; end
    case {'ls', 'jump_start', 'late_onset'}
        nq = 2; thresholds = -0.5:1:1.5;
    case 'moving_flies'
        nq = 5; thresholds = -0.5:1:5.5;
    otherwise
        nq = 4;
        thresholds = quantile(T.(param), linspace(0, 1, nq + 1));
end
end

function cmap = get_var_cmap(col, varname)
if isfield(col.vars, varname)
    cmap = col.vars.(varname);
else
    cmap = col.vars.default;
end
end

function labels = level_labels(param, nq, T)
switch param
    case {'nloom', 'ln'}
        labels = {sprintf('%d–5', max(1, min(T.nloom))), '6–10', '11–15', '16–20'};
    case 'ls'
        labels = {'Slow', 'Fast'};
    case 'late_onset'
        labels = {'Early', 'Late'};
    case 'moving_flies'
        labels = arrayfun(@(x) sprintf('%d moving', x), 0:nq-1, 'UniformOutput', false);
    otherwise
        labels = arrayfun(@(x) ['Q' num2str(x)], 1:nq, 'UniformOutput', false);
end
end
