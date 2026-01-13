function fh = plot_estimates(varargin)

opt = inputParser;
addParameter(opt, 'results', []);
addParameter(opt, 'export', false);
addParameter(opt, 'ylimits', [-1, 4]);
addParameter(opt, 'base', []);
addParameter(opt, 'paths', []);
addParameter(opt, 'marker', 'o');
addParameter(opt, 'label', []);
addParameter(opt, 'colors', 'flat');

parse(opt, varargin{:});
export = opt.Results.export;
paths = opt.Results.paths;
ylimits = opt.Results.ylimits;
figure_handle = opt.Results.base;
marker = opt.Results.marker;
label = opt.Results.label;
colors =  opt.Results.colors;
results = opt.Results.results;
est_means = results.estimates_mean(:, ~ismissing(results.estimates_mean));
est_std = results.estimates_std(:, ~ismissing(results.estimates_mean));


if isempty(figure_handle)

    quantiles = 0;
    col = cmapper([], quantiles);

    fh = figure('color','w', 'Position', [100 100 800 400]);

    hold on
    [suffixes, prefixes] = extract_dep(est_means);

    % Replace 'intercept' with '0' in the suffixes
    suffixes_replaced = suffixes;
    suffixes_replaced(strcmp(suffixes, 'intercept')) = {'0'};

    % Combine using cellfun
    xx = 1:size(suffixes, 2);

    result = cellfun(@(pre, suf) ['$$\beta_{' pre '}^{' suf '}$$'], ...
        prefixes, suffixes_replaced, 'UniformOutput', false);
    c = arrayfun(@num2str, xx, 'UniformOutput', false);
    result = cellfun(@(suf) ['$$\beta^{' suf '}$$'], ...
        c, 'UniformOutput', false);

    for idx_param = 1:length(xx)
    fill([xx(idx_param) - 0.3, xx(idx_param) - 0.3, xx(idx_param) + 0.3, xx(idx_param) + 0.3], ...
        [ylimits, fliplr(ylimits)], '','FaceColor', col.vars.(suffixes{idx_param}), 'LineStyle', 'none', 'FaceAlpha', 0.3,'HandleVisibility','off');
    end

    xticklabels(result);
else
    figure(figure_handle)
end


est_means = table2array(est_means);
est_std = table2array(est_std);
xx = 1:size(est_means, 2);


if isempty(figure_handle)
    errbar(xx, est_means, est_std, 'color','k','Linewidth', 1,'HandleVisibility','off');
    scatter(xx, est_means, 130, 'Marker', marker, 'MarkerFaceColor', 'k', 'HandleVisibility','off', 'MarkerEdgeColor', 'none');
else
    sh = scatter(xx, est_means, 130, 'Marker', marker, 'MarkerFaceColor', 'none', 'DisplayName', label, 'LineWidth', 1.3, 'MarkerEdgeColor', colors);
    errbar(xx, est_means, est_std,'Linewidth', 1.1,'HandleVisibility','off', 'Color', sh.MarkerEdgeColor);
end

plot([xx(1) - 1,xx(end) + 1], [0 0], 'k--','HandleVisibility','off');

if isfield(results, 'ground_truth')
    ground_truth = table2array(results.ground_truth(1, ~ismissing(results.estimates_mean)));
    ground_truth(isnan(ground_truth)) = 0;
    scatter(xx, ground_truth, 200, 'k', 'Marker', 'diamond', 'LineWidth', 1)
end

xlim([xx(1) - 1,xx(end) + 1]); ylim(ylimits)
ax = gca;
apply_generic(ax)
xticks(xx);
set(ax.XAxis, 'TickLabelInterpreter', 'latex', 'FontSize', 24);


if export
    paths.fig = results.fig_path;
    exporter(fh, paths, 'estimates.pdf')
end


function [suffixes, prefixes] = extract_dep(results)

params = results.Properties.VariableNames;
parts_str = cellfun(@(s) split(s, '_'), params, 'UniformOutput', false);
suffixes = cellfun(@(parts) parts{end}, parts_str, 'UniformOutput', false);
prefixes = cellfun(@(concat) concat{1:end-1}, parts_str, 'UniformOutput', false);