
function overlay_fits(fh, ax, ax_inset, varargin)

opt = inputParser;

addParameter(opt, 'extra', []);
addParameter(opt, 'freezes', []);
addParameter(opt, 'results', []);
addParameter(opt, 'bin_size', 3);
addParameter(opt, 'conditions', false);
addParameter(opt, 'type', 'continuous');
addParameter(opt, 'color', 'col');
addParameter(opt, 'censored_inset', true);
addParameter(opt, 'gt', false);
addParameter(opt, 'export', false);
addParameter(opt, 'paths', []);

parse(opt, varargin{:});

extra = opt.Results.extra;
freezes = opt.Results.freezes;
results = opt.Results.results;
conditions = opt.Results.conditions;
export = opt.Results.export;
paths = opt.Results.paths;
censored_inset = opt.Results.censored_inset;
bin_size = opt.Results.bin_size;
gt_plot = opt.Results.gt;
type = opt.Results.type;
color = opt.Results.color;

est_params = table2array(results.estimates_mean(:, ~ismissing(results.estimates_mean)));
freezes = importdata(fullfile(results.bouts_path, 'surrogate.mat'));

i = 0;
for idx_sm = 1:3
    for idx_ls = 1:2
        for idx_fs = 1:2

            i = i + 1;
            axes(ax(i));

            [freezes_quant, ~] = quantilizer(freezes, 'idx_quanti', struct('sm', idx_sm, 'ls', idx_ls, 'fs', idx_fs));

            [~, f, fd] =  nll_fly_ddm_newer(est_params, freezes_quant, results.points, strcat('model_', results.fitted_model), 'iid', 'p', extra);

            plot(fd, f, 'LineWidth', 1.25, 'Color', 'k', 'LineStyle', '--')

            axes(ax_inset(i));
            plot(results.points.censoring, f(end), 'o', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerEdgeColor', 'k');

        end
    end
end

if export
    paths.fig = results.fig_path;
    exporter(fh, paths, 'fits.pdf')
end