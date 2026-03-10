
% Load the table with the computed log likelihoods 
model_list = {'dddim2'};
run_list   = {'run01_010326_1800'};
version = {'_v8'};

paths_analysis = path_generator('folder', fullfile('momentary_integration','alternative_test', model_list{1}, strcat(run_list{1}, version{1})));
paths_analysis.fig = fullfile(paths_analysis.fig);

col = cmapper('', 2);

ll = struct2table(importdata(fullfile(paths_analysis.results, 'll_table.mat')));
data = alternative_test_collect_data(model_list, run_list);

fh = figure('color','w','Position',[100, 100, 1200, 500]);
tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact')
nexttile

total_ll_bar(ll, col)

nexttile(2, [1 2])
hold on
ylim([-2 2])
xlim([-50 height(ll) + 50])
ax = gca;

[sorted_deltall, idx_deltall_tv] = sort(ll.tv - ll.st);
idx_crossing = find(sorted_deltall > 0);

fill([ax.XLim(1) ax.XLim(1) idx_crossing(1) idx_crossing(1)], [ax.YLim(1) 0 0 ax.YLim(1)], hex2rgb(col.stationary_sm), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'LineWidth', 2)
fill([idx_crossing(1) idx_crossing(1) ax.XLim(2) ax.XLim(2)], [ax.YLim(2) 0 0 ax.YLim(2)], hex2rgb(col.timevarying_sm), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'LineWidth', 2)

xline(idx_crossing(1), 'Label', sprintf('%.2f%%', (sum(sorted_deltall > 0)./length(sorted_deltall)) * 100), 'FontSize', 20, 'LabelOrientation','horizontal')

bar(sorted_deltall, 'FaceColor', 'k', 'EdgeColor', 'none' )
xlabel('Sorted Freezes')
xticks([])
ylabel('$\log \!\big(\mathcal{L}_{\mathrm{tv}} / \mathcal{L}_{\mathrm{st}}\big)$', ...
    'Interpreter','latex')
apply_generic(ax)

nexttile(4)
hold on
bin_size = 0.020;
histogram(sorted_deltall(sorted_deltall>0), 0:bin_size:3, 'Orientation', 'horizontal', 'EdgeColor', 'none', 'FaceColor', hex2rgb(col.timevarying_sm))
histogram(sorted_deltall(sorted_deltall<0), -3:bin_size:0, 'Orientation', 'horizontal', 'EdgeColor', 'none', 'FaceColor', hex2rgb(col.stationary_sm))

xlabel('Count')
apply_generic(gca, 'no_y', true, 'ylim', [-2 2])

exporter(fh, paths_analysis, 'compare_total_ll.pdf')

function total_ll_bar(lls_output, col)

bh = bar([sum(lls_output.tv), sum(lls_output.st)], 'FaceColor', 'flat', 'EdgeColor', 'flat', 'LineWidth', 2);
bh.CData(1,:) = hex2rgb(col.timevarying_sm);
bh.CData(2,:) = hex2rgb(col.stationary_sm);

xticklabels({'time-varying', 'constant'})
ylabel('Total log($\mathcal{L}$)', 'Interpreter', 'latex')
ax = gca;
apply_generic(ax)
xtickangle(25)
ax.YAxis.Direction = 'reverse'; %ax.YTickLabel(1) = {'worse'}; ax.YTickLabel(end) = {'better'};
text(2, ax.YLim(1) + 500, sprintf('$\\Delta log(\\mathcal{L}): %.2f$', sum(lls_output.tv) - sum(lls_output.st)), 'FontSize', 20, 'Interpreter', 'latex', 'HorizontalAlignment', 'center')

end
