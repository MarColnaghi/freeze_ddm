function plot_aligned_offset(offset_signal, al, sort_order, bouts_proc, colors, paths, export_name, show_offset_xline)
% PLOT_ALIGNED_OFFSET  Two-panel figure used for both the social-motion and the
% min-distance versions (offset_signal is the only thing that differs).
%   Left panel : peak-aligned social motion split into before/during/after the
%                freeze break (always al.aligned_peak_valid).
%   Right panel: the offset-aligned OFFSET_SIGNAL split by pre-break duration,
%                with a zoom inset.
%
%   export_name        '' to skip export, else base filename (no extension)
%   show_offset_xline   draw the dashed x=0 line on the right panel

if nargin < 8 || isempty(show_offset_xline), show_offset_xline = true; end
if nargin < 7, export_name = ''; end

aligned_peak_valid     = al.aligned_peak_valid;
x_axis_peak            = al.x_axis_peak;
break_times_peak_valid = al.break_times_peak_valid;
valid_sorted           = al.valid_sorted;

delay_after = 60;

fh = figure('color','w','Position',[100 100 750 400]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact')

% ---- LEFT: peak-aligned SM, before/during/after the break ----
nexttile
hold on
window = 15;

before = break_times_peak_valid < -window;
after  = break_times_peak_valid > window;
during = ~(after | before);

conds = {before, during, after};

for i = 1:3
    idxc = conds{i};
    data = aligned_peak_valid(idxc, :);

    mu  = mean(data, 1, 'omitnan');
    sem = std(data, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data), 1));
    m = ~isnan(sem);

    if i == 1   % shaded error only for "before"
        fill_between(x_axis_peak(m), mu(m) + sem(m), mu(m) - sem(m), [], ...
            'FaceColor', colors(i, :), 'EdgeColor', 'none', 'FaceAlpha', 0.25);
    end

    plot(x_axis_peak, mu, 'Color', colors(i, :), 'LineWidth', 2);
end

apply_generic(gca, 'xlim', [-360 120]./2, 'xticks', -180:60:180, 'ylim', [0 2]);
xticklabels(-3:1:3)
ylabel('Social Motion')
xlabel('Time (Aligned to Peak)')

% ---- RIGHT: offset-aligned signal, split by pre-break duration ----
nexttile
hold on
if show_offset_xline
    xline(0, 'k-.')
end

signal_valid = offset_signal(sort_order(valid_sorted), :);
durs = bouts_proc.durations(sort_order(valid_sorted));

before_early = before & durs >= 30 & durs < 60;
before_mid   = before & durs >= 60 & durs < 180;
before_late  = before & durs >= 180;

sm_sets = {
    signal_valid(before_early,:)
    signal_valid(before_mid,:)
    signal_valid(before_late,:)
    signal_valid(during,:)
    signal_valid(after,:)
    };

colors_before = cbrewer2('Blues', 4);
colors_loop = [colors_before(2:end,:); colors(2:3, :)];

for i = 1:5
    data = sm_sets{i};

    mu  = mean(data, 1, 'omitnan');
    sem = std(data, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data), 1));
    m = ~isnan(sem);

    t = -numel(mu) + 1 + delay_after : delay_after;

    fill_between(t(m), mu(m) + sem(m), mu(m) - sem(m), [], ...
        'FaceColor', colors_loop(i, :), 'EdgeColor', 'none', 'FaceAlpha', 0.25);

    plot(t, mu, 'Color', colors_loop(i, :), 'LineWidth', 2);
end

zoomY = [0.3 .8];
zoomX = [-120 delay_after];

fill([zoomX fliplr(zoomX)], [zoomY(1) zoomY(1) zoomY(2) zoomY(2)], [0 0 0], 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 2)
apply_generic(gca, 'xlim', [-300 delay_after], 'xticks', -360:120:0, 'ylim', [0 2])
xticklabels(-6:2:0)
xlabel('Time (Aligned to Offset)')

% ---- zoom inset on the right panel ----
ax2 = gca;
pos = ax2.Position;
inset_pos = [pos(1) + 0.2*pos(3), pos(2) + 0.55*pos(4), 0.5*pos(3), 0.35*pos(4)];

ax_inset = axes('Position', inset_pos);
hold(ax_inset, 'on')
box(ax_inset, 'on')

for i = 1:3
    data = sm_sets{i};
    if isempty(data)
        continue
    end

    mu  = mean(data, 1, 'omitnan');
    sem = std(data, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(data), 1));
    m = ~isnan(sem);
    t = -numel(mu) + 1 + delay_after : delay_after;

    fill_between(t(m), mu(m) + sem(m), mu(m) - sem(m), [], ...
        'FaceColor', colors_loop(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.25);

    plot(t, mu, 'Color', colors_loop(i,:), 'LineWidth', 1.5);

    apply_generic(ax_inset, 'xlim', zoomX, 'xticks', -360:120:0, 'ylim', zoomY, 'yticks', [0.3 .8])
    xticklabels(-6:2:0)
end

if ~isempty(export_name)
    exporter(fh, paths, [export_name '.pdf'])
    exporter(fh, paths, [export_name '.png'])
end
end
