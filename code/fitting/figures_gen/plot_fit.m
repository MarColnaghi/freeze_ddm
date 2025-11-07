function plot_fit(varargin)

opt = inputParser;
addParameter(opt, 'extra', []);
addParameter(opt, 'freezes', []);
addParameter(opt, 'results', []);
addParameter(opt, 'conditions', false);
addParameter(opt, 'export', false);
addParameter(opt, 'paths', []);

parse(opt, varargin{:});

extra = opt.Results.extra;
freezes = opt.Results.freezes;
results = opt.Results.results;
conditions = opt.Results.conditions;
export = opt.Results.export;
paths = opt.Results.paths;

if isempty(opt.Results.freezes)
     load(fullfile(opt.Results.results.bouts_path, 'surrogate.mat'));
     freezes = surrogate;

end

est_params = table2array(results.estimates_mean(:, ~ismissing(results.estimates_mean)));
results.fitted_model = sprintf('model_%s', results.fitted_model);

if ~conditions

    fh = figure('Position', [100 100 800 400], 'Color', 'w');
    [~, f, fd] = nll_fly_ddm_newer(est_params, freezes, results.points, results.fitted_model, 'iid', 'p', extra);
    hold on
    histogram(freezes.durations_s, -1/120:1/60:12, 'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none')
    plot(fd,f, 'k--', 'LineWidth', 1.5)
    apply_generic(gca)
    xlabel('Freeze Duration (s)')
    ylabel('pdf')

else
    fh = figure('Position', [100 100 1400 900], 'Color', 'w');
    tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'loose')
    i = 0;
    for idx_sm = 1:3
        for idx_ls = 1:2
            for idx_fs = 1:2

                i = i + 1;
                nexttile
                hold on

                freezes_quant = quantilizer(freezes, 'idx_quanti', struct('sm', idx_sm, 'ls', idx_ls, 'fs', idx_fs));
                [~, f, fd] = nll_fly_ddm_newer(est_params, freezes_quant, results.points, results.fitted_model, 'iid', 'p', extra);

                histogram(freezes_quant.durations_s, -1/120:1/15:12, 'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none')
                plot(fd,f, 'k--', 'LineWidth', 1.5)
                ax(i) = gca;
                apply_generic(gca)
                xlim([-0.1 10.1])
                % ylim([-0.03 1.53])
               
            end
        end
    end
    linkaxes(ax)

end
if export
    exporter(fh, paths, 'fits.pdf')
end
