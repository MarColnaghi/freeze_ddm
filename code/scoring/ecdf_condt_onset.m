%% Script Name:
% ecdf_condt_onset
% last update: 09.03

function ecdf_condt_onset(bouts, thresholds)

% Load Colors
extra.quantiles = 5;
col = cmapper([], extra.quantiles);

if nargin == 0
    %clear all; close all; clc

    disp('no input provided, will load the dataset in the dataset folder')
    genotype = 1;

    % Load Paths & Files
    paths = path_generator('validation_scoring/durs_onsets');
    load(fullfile(paths.dataset, 'bouts.mat'));
    bouts = bouts(bouts.genotype == genotype, :);

    % Load Colors
    extra.quantiles = 61;
    col = cmapper([], extra);

    zoom = true;


    % Define thresholds
    thresholds = define_thresholds;
    time_window = thresholds.time_window;
    bouts = bouts_formatting(bouts, thresholds);


else
    disp('input provided')
    time_window = thresholds;
end

% Plot Histogram 2D of Freeze Onset and Freeze Duration
str = 'immobility';
type = 'loom';

% Condition Flags
clamp = false;
bouts_with_loom = 'all';
zoom =  false;

% Number of onsets
n_onsets = 240;
start_position = - 10;
fh = figure('color','w','Position',[100, 100, 1200, 600]);
tl = tiledlayout(1,2, "TileSpacing","compact","Padding","compact");

i = 0;
for idx_slooms = [unique(bouts.sloom)]'
    i = i + 1;
    
    nexttile
    if idx_slooms == 25
        str = 'le_window_sl';
    elseif idx_slooms == 50
        str = 'le_window_fl';
    end

    bouts_temp = bouts(bouts.sloom == idx_slooms, :);

    if clamp == true
        bouts_temp.durations(bouts_temp.durations + bouts_temp.onsets_loomwin >= time_window(2)) = time_window(2) - bouts_temp.onsets_loomwin(bouts_temp.durations + bouts_temp.onsets_loomwin >= time_window(2));
    end

    if islogical(bouts_with_loom)
        bouts_temp = bouts_temp(bouts_temp.bout_with_loom == bouts_with_loom, :);
    end

    if strcmp(str, 'immobility')
        bouts_temp = bouts_temp(bouts_temp.type == 1, :);

    elseif strcmp(str, 'mobility')
        bouts_temp = bouts_temp(bouts_temp.type == 0, :);
    end

    % Create Figure
    hold on
    cmap = colorcet('C06', 'N', n_onsets + abs(start_position) + 1);

    for idx_onsets = start_position:1:n_onsets
        [f, x] = ecdf(bouts_temp.durations(bouts_temp.onsets_loomwin == idx_onsets)); %, 'censoring', bouts_temp.bout_with_loom(bouts_temp.onsets_loomwin == idx_onsets));

        if idx_onsets >= thresholds.(str)(1) && idx_onsets <= thresholds.(str)(2)
            plot(x, f, 'Color', cmap(idx_onsets + abs(start_position) + 1, :) , 'LineWidth', 1)
        else
            plot(x, f, 'Color', 'k' , 'LineWidth', 1)
        end

    end
    ax(i) = gca;
    ylim([0 1])
    xlim([0 1200])
    apply_generic(gca)
    xlabel('Immobility Duration (frames)')
    ylabel('ecdf')

    cb = colorbar('LineWidth', 2, 'Location', 'northoutside');
    set(gca, 'CLim', [0, size(cmap,1)]);
    cmap(1:thresholds.(str)(1)+ abs(start_position) + 1, :) = repmat(zeros(1,3), size(cmap(1:thresholds.(str)(1)+ abs(start_position) + 1, :), 1), 1);
    cmap(thresholds.(str)(2) + abs(start_position) + 1:end, :) = repmat(zeros(1,3), size(cmap(thresholds.(str)(2) + abs(start_position) + 1:end, :), 1), 1);
    colormap(ax(i), cmap);

    cb.Label.String = 'Immobility Onset';
    cb.FontSize = 20;
    cb.Ticks = [0 30 60 120 240] + abs(start_position);
    cb.TickLabels = {'on', 'off', 60, 120, 240};

end
