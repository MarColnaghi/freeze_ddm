
function overlay_fits_for_cosyne(fh, ax, ax_inset, varargin)

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
% est_params = [  -1.0101       3.0639        1.2792      0.029802       -0.017144        0.24792         -0.52892          2.8024               3               0.0001    ];

freezes = importdata(fullfile(results.bouts_path, 'freeze.mat'));

i = 0;
for idx_sm = 1:3
    for idx_ls = 0:1
        for idx_fs = 1:2

            ec = extra;

            i = i + 1;
            axes(ax(i));

            [freezes_quant, mask] = quantilizer_v2(freezes, 'indexed_quantile', struct('sm', idx_sm, 'ls', idx_ls, 'fs', idx_fs));

            if isfield(extra, 'soc_mot_array')
                ec.soc_mot_array = extra.soc_mot_array(mask, :);
            end

            [~, f, fd] =  nll_fly_ddm_newer(est_params, freezes_quant, results.points, strcat('model_', results.fitted_model), 'iid', 'p', extra);

            plot(fd, f, 'LineWidth', 1.9, 'Color', 'k', 'LineStyle', '--')

            if censored_inset
                axes(ax_inset(i));
                plot(results.points.censoring, f(end), 'o', 'LineWidth', 1, 'MarkerSize', 5, 'MarkerEdgeColor', 'k');
            end

        end
    end
end

if export
    paths.fig = results.fig_path;
    exporter(fh, paths, sprintf('fits_%s.pdf', color))
end