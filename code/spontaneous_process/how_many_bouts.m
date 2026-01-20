% This script shows how many bouts are there per fly, and a summary of the
% distribution of counts at the population level
% 

id_code = 'imm2_mob2_pc4';
paths_out = path_generator('folder', '/spontaneous_process', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));

col = cmapper('', 5);
fh = figure('Position', [100 100 800 410], 'Color', 'w');
t = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'loose');
bouts_spontaneous = data_parser_new(bouts, 'period', 'bsl', 'window', 'all', 'type', 'immobility');

[counts_x_fly] = histcounts(bouts_spontaneous.fly, 0.5:1:max(bouts_spontaneous.fly) + 0.5);
[~, ~, fly_idx] = unique(bouts_spontaneous.fly, 'stable');
fly_moving = accumarray(fly_idx, bouts_spontaneous.moving_flies, [], @mode);
fly_colors = col.vars.moving_flies(fly_moving + 1, :);

[sorted_counts, sorted_i] = sort(counts_x_fly);
sorted_colors = fly_colors(sorted_i, :);

nexttile
bar(sort(counts_x_fly), 'EdgeColor', 'none', 'FaceColor', 'flat', 'CData', sorted_colors)
apply_generic(gca, 'tick_length', 0.01)
xlabel('Flies')
ylabel('Counts')

nexttile
histogram(counts_x_fly, 50, 'EdgeColor', 'none')
apply_generic(gca, 'tick_length', 0.01)
xlabel('Number of Bouts')
ylabel('Counts')

exporter(fh, paths_out, 'contaminant_howmany_xfly.pdf')