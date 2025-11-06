%% Script Name:
% hist2d_durs_onsets
% last update: 04.03

function hist2d_durs_onsets(bouts, thresholds, type, str, plot_style, zoom_flag, paths, export)

% Load Colors
extra.quantiles = 5;
col = cmapper([], extra);

if nargin == 0
    %clear all; close all; clc

    disp('no input provided, will load the dataset in the dataset folder')

    genotype = 1;

    % Load Paths & Files
    paths = path_generator('validation_scoring/durs_onsets');
    load(fullfile(paths.dataset, 'bouts.mat'));
    bouts = bouts(bouts.genotype == genotype, :);

    zoom_in = 100;

    % Condition Flags
    clamp = false;
    bouts_with_loom = 'all';
    zoom =  false;

    % Define thresholds
    thresholds = define_thresholds;
    time_window = thresholds.time_window;
    bouts = bouts_formatting(bouts, thresholds);

    % Plot Histogram 2D of Freeze Onset and Freeze Duration
    str = 'immobility';
    type = 'loom';
    plot_style = 'tile';
    export = false;

else
    disp('input provided')
    time_window = thresholds.time_window;
    
    clamp = false;
    bouts_with_loom = false;
    zoom =  zoom_flag;
    limits = [0 100];
end

switch plot_style
    case 'bar3'
        bin_size = 1;
    case 'tile'
        bin_size = 2;
end

i = 0; 
for idx_slooms = [unique(bouts.sloom)]'
    i = i + 1;
    bouts_temp = bouts;
    bouts_temp = bouts(bouts.sloom == idx_slooms, :);
    % bouts_temp = bouts(bouts.loom_durs == 29, :);

    if clamp == true
        bouts_temp.durations(bouts_temp.durations + bouts_temp.onsets_loomwin >= time_window(2)) = time_window(2) - bouts_temp.onsets_loomwin(bouts_temp.durations + bouts_temp.onsets_loomwin >= time_window(2));
    end

    if islogical(bouts_with_loom)
        bouts_temp = bouts_temp(bouts_temp.bout_with_loom == bouts_with_loom, :);
    end

    if strcmp(type, 'loom')
        bouts_temp = bouts_temp(bouts_temp.period == 1, :);

    elseif strcmp(type, 'bsl')
        bouts_temp = bouts_temp(bouts_temp.period == 0, :);
    end

    if strcmp(str, 'immobility')
        bouts_temp = bouts_temp(bouts_temp.type == 1, :);

    elseif strcmp(str, 'mobility')
        bouts_temp = bouts_temp(bouts_temp.type == 0, :);
    end

    % Create Figure
    fh = figure('color','w','Position',[100, 100, 800, 800]);
    hold on

    h2h = histogram2(bouts_temp.durations, bouts_temp.onsets_loomwin, -0.5:bin_size:max(bouts_temp.durations) + 0.5, -time_window(1) - 0.5:bin_size:time_window(2) + 0.5, 'DisplayStyle', plot_style, 'FaceColor', 'flat');
    % If you want marginals
    % h2h = histogram2(durs, 604 * ones(length(durs), 1), -0.5:1:601.5, 603:605, 'DisplayStyle', 'bar3', 'FaceColor', 'k', 'Normalization', 'pdf', 'EdgeColor','none');
    % h2h = histogram2(-1 * ones(length(onsets), 1), onsets, -2:0,  -0.5:1:601.5, 'DisplayStyle', 'bar3', 'FaceColor', 'k', 'Normalization', 'pdf', 'EdgeColor','none');

    % Set Axis Labels
    xlabel('Bout Duration (frames)');
    ylabel('Bout Latency (frames)');
    zlabel('Count');

    % Set Axis Limits
    ylim([-time_window(1) - 10 time_window(2) + 25])
    xlim([-10 710]);
    zlim([-0.2, max(h2h.Values, [], 'all') + 10]);

    if zoom
        ylim([-10 61])
        xlim([-1 81]);
        zlim([-0.2, max(h2h.Values, [], 'all') + 10]);
    end

    % Set Colormap
    colorcet('L08');
    clim([0 max(h2h.Values, [], 'all')])

    % Apply Generic Changes
    ax = gca;
    apply_generic(ax, 20)
    set(ax,'TickLength',[0.02, 0.02])
    grid on

    switch plot_style
        case 'bar3'
            thresholds.sp_window = [];
            ax.XLabel.Rotation = 300;
            ax.XLabel.Position(1:2) = [250 -225];
            ax.YLabel.Rotation = 7;
            ax.YLabel.Position(1:2) = [675 210];
            if ~zoom
                view(80, 50)
            else
                view(80, 30)

            end
            % Set Colormap
            colorcet('L08');
            clim([0 round(max(h2h.Values, [], 'all'), -1)])

        case 'tile'
            bin_size = 2;
            view(90, 90)
            % Set Colormap
            colorcet('L02');
            clim([0 round(max(h2h.Values, [], 'all'), -1)])

    end

    % Set Axis Ticks
    if ~zoom
        xticks(0:100:600)
    else
        xticks([0 60])
    end
    ticks = sort([-time_window(1) 0 30 100 200 300 400 500 time_window(2)]);
    yticks(ticks)
    ytlh = yticklabels; ytlh{2} = 'onset'; ytlh{3} = 'offset'; yticklabels(ytlh); ytickangle(360 - 60);

    % Create Plane of Loom Presentation
    x = meshgrid(-10:905:900);
    y = meshgrid(0:30:30)';
    z = [-0.1 -0.1; -0.1 -0.1];

    zmin = ax.ZLim(1);
    zmax = ax.ZLim(2);
    zrange = zmax - zmin;

    z_label_pos = zmax + 0.05 * zrange;  % 5% above the top of the plot

    switch type
        case 'loom'
            if idx_slooms == 25
                lp = surf(x, y, z, 'FaceColor', col.vars.sloom(3,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'LineWidth', 1.5);
                surf(z, y, x, 'FaceColor', col.vars.sloom(3,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'LineWidth', 1.5);
                text(0, 15, z_label_pos, {'Slow', 'Loom'}, 'FontSize', 18, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Rotation', ax.YLabel.Rotation)
            elseif idx_slooms == 50
                lp = surf(x, y, z, 'FaceColor', col.vars.sloom(5,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'LineWidth', 1.5);
                surf(z, y, x, 'FaceColor', col.vars.sloom(5,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'LineWidth', 1.5);
                text(0, 15, z_label_pos, {'Fast', 'Loom'}, 'FontSize', 18, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Rotation', ax.YLabel.Rotation)

            end
        case 'bsl'
            text(0, 15, z_label_pos, {'Surrogate', 'Loom'}, 'FontSize', 18, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Rotation', ax.YLabel.Rotation)

    end
    
   % if zoom

        if idx_slooms == 25
            x = [0 0; 900 900];
            y_ons = meshgrid(linspace(thresholds.le_window_sl(1), thresholds.le_window_sl(1), 2))';
            y_off = meshgrid(linspace(thresholds.le_window_sl(2), thresholds.le_window_sl(2), 2))';

            z = [0 900; 0 900];
            surf(x, y_ons, z, 'FaceColor', 'k', 'FaceAlpha', 0.05, 'EdgeColor', 'k', 'LineWidth', 1.5);
            surf(x, y_off, z, 'FaceColor', 'k', 'FaceAlpha', 0.05, 'EdgeColor', 'k', 'LineWidth', 1.5);

        elseif idx_slooms == 50
            x = [0 0; 900 900];
            y_ons = meshgrid(linspace(thresholds.le_window_fl(1), thresholds.le_window_fl(1), 2))';
            y_off = meshgrid(linspace(thresholds.le_window_fl(2), thresholds.le_window_fl(2), 2))';

            z = [0 900; 0 900];
            surf(x, y_ons, z, 'FaceColor', 'k', 'FaceAlpha', 0.05, 'EdgeColor', 'k', 'LineWidth', 1.5);
            surf(x, y_off, z, 'FaceColor', 'k', 'FaceAlpha', 0.05, 'EdgeColor', 'k', 'LineWidth', 1.5);

        end
    %end

    if ~isempty(thresholds.sp_window)
        lp.Visible = 'off';
        for idx_window = 1:size(thresholds.sp_window, 1)
            if thresholds.sp_window(idx_window, 1) == 45
                col.areas(end - idx_window, :) = [87 167 115]/255;
            end
            surf(meshgrid(-5:sum(time_window):sum(time_window)), meshgrid(thresholds.sp_window(idx_window, 1):30:thresholds.sp_window(idx_window, 2))',...
                [-0.01 -0.01; -0.01 -0.01], 'FaceColor', col.areas(end - idx_window, :), 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'LineWidth', 1);
        end
    end

    % Export Figure
   
    if export
        figure_title = sprintf('%s_%s_%s_sloom%d_z%d', plot_style, str, type, idx_slooms, zoom_flag);

        exporter(ax, paths, [figure_title, '.png'])
        exporter(ax, paths, [figure_title, '.pdf'])
    end
end