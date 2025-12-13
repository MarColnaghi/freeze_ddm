function [fh] = base_for_estimates(varargin)

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