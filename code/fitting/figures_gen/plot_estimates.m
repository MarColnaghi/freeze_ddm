function plot_estimates(varargin)

opt = inputParser;
addParameter(opt, 'results', []);
addParameter(opt, 'export', false);
addParameter(opt, 'ylimits', [-1, 4]);

addParameter(opt, 'paths', []);

parse(opt, varargin{:});
export = opt.Results.export;
paths = opt.Results.paths;
ylimits = opt.Results.ylimits;

quantiles = 0;
col = cmapper([], quantiles);

results = opt.Results.results;
est_means = results.estimates_mean(:, ~ismissing(results.estimates_mean));
est_std = results.estimates_std(:, ~ismissing(results.estimates_mean));

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

est_means = table2array(est_means);
est_std = table2array(est_std);

errbar(xx, est_means, est_std, 'color','k','Linewidth', 2,'HandleVisibility','off');
scatter(xx, est_means, 100, 'filled', 'MarkerFaceColor', 'k', 'HandleVisibility','off');
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
xticklabels(result); set(ax.XAxis, 'TickLabelInterpreter', 'latex', 'FontSize', 24);

if export
    paths.fig = results.fig_path;
    exporter(fh, paths, 'estimates.pdf')
end


function [suffixes, prefixes] = extract_dep(results)

params = results.Properties.VariableNames;
parts_str = cellfun(@(s) split(s, '_'), params, 'UniformOutput', false);
suffixes = cellfun(@(parts) parts{2}, parts_str, 'UniformOutput', false);
prefixes = cellfun(@(parts) parts{1}, parts_str, 'UniformOutput', false);