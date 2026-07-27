% This script shows how many bouts are there per fly, and a summary of the
% distribution of counts at the population level


id_code = 'imm2_mob2_pc4';
paths_out = path_generator('folder', 'spontaneous_process', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));

col = cmapper('', 5);
fh = figure('Position', [100 100 800 410], 'Color', 'w');
t = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'loose');
bouts_spontaneous = data_parser_new(bouts, 'period', 'bsl', 'type', 'immobility', 'window', 'all');

[~, ~, fly_idx] = unique(bouts_spontaneous.fly, 'stable');
% Count bouts per fly using the same basis as the per-fly stats below, so
% indexing stays aligned even when fly ids are non-contiguous.
counts_x_fly = accumarray(fly_idx, 1)';

select_method = 'median';

% Median duration per fly
dur_per_fly.median = accumarray( ...
    fly_idx, ...
    bouts_spontaneous.durations_s, ...
    [], ...
    @median);

dur_per_fly.mean = accumarray( ...
    fly_idx, ...
    bouts_spontaneous.durations_s, ...
    [], ...
    @mean);

dur_per_fly.sum = accumarray( ...
    fly_idx, ...
    bouts_spontaneous.durations_s, ...
    [], ...
    @sum);


fly_moving = accumarray(fly_idx, bouts_spontaneous.moving_flies, [], @mode);
fly_colors = col.vars.moving_flies(fly_moving + 1, :);

[sorted_counts, sorted_i] = sort(counts_x_fly);
sorted_colors = fly_colors(sorted_i, :);
sorted_medians = dur_per_fly.(select_method)(sorted_i);

nexttile
%yyaxis left

bar(sort(counts_x_fly), 'EdgeColor', 'none', 'FaceColor', 'flat', 'CData', sorted_colors)
xlabel('Flies')
ylabel('Counts')
ax = gca;
apply_generic(ax)
hold on
yyaxis right
plot(1:length(sorted_medians), repmat(0.1, length(sorted_medians), 1), 'r-.')
plot(sorted_medians, 'r')
ylabel('Median Duration', 'FontSize', 20)

nexttile
histogram(counts_x_fly, 50, 'EdgeColor', 'none')
apply_generic(gca, 'tick_length', 0.01)
xlabel('Number of Bouts')
ylabel('Counts')

exporter(fh, paths_out, 'contaminant_howmany_xfly.pdf')

