function compare_elbos(varargin)

opt = inputParser;
addParameter(opt, 'model', []);
addParameter(opt, 'run', []);
addParameter(opt, 'export', false);
addParameter(opt, 'ylimits', [-1, 4]);

addParameter(opt, 'paths', []);

parse(opt, varargin{:});
model_list = opt.Results.model;
run_list = opt.Results.run;
export = opt.Results.export;

quantiles = 0; i = 0;
col = cmapper([], quantiles);
col.Set3 = cbrewer2('Set3', 5);

paths = path_generator('folder', 'fitting_freezes/le');

fh = figure('Position', [100 100 400 700], 'Color', 'w');
t = tiledlayout(1, 1, 'TileSpacing', 'loose', 'Padding', 'loose');
nexttile
ax = gca;
elbos = nan(length(model_list), 2);

for idx_models = 1:length(model_list)
    
    model = model_list{idx_models}; run = run_list{idx_models};
    model_func = str2func(sprintf('model_%s', model)); mod = model_func();
    model_results = importdata(fullfile(paths.results, model, run, sprintf('fit_results_%s.mat', model)));
    y = importdata(fullfile(paths.results, model, run, 'surrogate.mat'));
    param_table = model_results.estimates_mean;

    elbos(idx_models, :) = model_results.elbo;
    
end

bar(elbos(:,1), 'EdgeColor', 'none')
apply_generic(ax, 'ylim', [-12500 -12000])
ax.YAxis.Direction = 'reverse';
xticklabels({'Fixed Rate', 'Condition Dependent'})
paths.fig = fullfile(paths.fig);

if export
    exporter(fh, paths, 'elbo_comparison.pdf')
end

