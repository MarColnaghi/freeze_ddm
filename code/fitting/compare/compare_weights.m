function compare_weights(varargin)

opt = inputParser;
addParameter(opt, 'model', []);
addParameter(opt, 'run', []);
addParameter(opt, 'export', false);
addParameter(opt, 'ylimits', [-1, 4]);

addParameter(opt, 'paths', []);

parse(opt, varargin{:});
model_list = opt.Results.model;
run_list = opt.Results.run;

quantiles = 0; i = 0;
col = cmapper([], quantiles);
col.Set3 = cbrewer2('Set3', 5);

paths = path_generator('folder', 'fitting_freezes/le');

fh = figure('Position', [100 100 1400 900], 'Color', 'w');
t = tiledlayout(3, 4, 'TileSpacing', 'loose', 'Padding', 'loose');

for idx_tiles = 1:12
    i = i + 1;
    nt =  nexttile(t);
    ax(i) = gca;
end

for idx_models = 1:length(model_list)
    
    i = 0;

    model = model_list{idx_models}; run = run_list{idx_models};
    model_func = str2func(sprintf('model_%s', model)); mod = model_func();
    model_results = importdata(fullfile(paths.results, model, run, sprintf('fit_results_%s.mat', model)));
    y = importdata(fullfile(paths.results, model, run, 'surrogate.mat'));
    param_table = model_results.estimates_mean;

    for idx_sm = 1:3
        for idx_ls = 1:2
            for idx_fs = 1:2

                i = i + 1;
                axes(ax(i))
                hold on
                [y_quant, mask] = quantilizer(y, 'idx_quanti', struct('sm', idx_sm, 'ls', idx_ls, 'fs', idx_fs));
                out = evaluate_model(mod, param_table, y_quant);

                cell_array =     out.Properties.VariableNames;
                pmix_cells = cell_array(contains(cell_array, 'pmix'));
                pmix = table;

                height(y_quant)
                if length(pmix_cells) <1
                    pmix.s = out.pmix;
                    pmix.l = 1 - out.pmix;

                else
                    if ~any(strcmp(pmix_cells, 'pmix_contaminant'))

                        pmix.c = out.pmix_sum .* out.pmix_ratio;
                        pmix.s = out.pmix_sum .* (1 - out.pmix_ratio);
                        pmix.l = 1 - pmix.c - pmix.s;
                    
                    elseif any(strcmp(pmix_cells, 'pmix_contaminant'))

                        pmix.c = out.pmix_contaminant;
                        pmix_remaining = 1 - pmix.c;
                        pmix.s = pmix_remaining .* out.pmix_ratio;
                        pmix.l = pmix_remaining .* (1 - out.pmix_ratio);

                    end

                end
                
                base = 1:size(pmix, 2);
                width = 0.15;
                offsets = linspace(-width, +width, length(model_list));
                xx = base + offsets(idx_models);
                %
                hand = plot(xx, pmix.Variables, 'o', 'MarkerFaceColor', col.Set3(idx_models + 2, :),...
                    'MarkerEdgeColor', 'none', 'MarkerSize', 1);

%                 hand = scatter(xx, pmix.Variables, 5, 'o', 'MarkerFaceColor', col.Set3(idx_models + 2, :),...
%                     'MarkerEdgeColor', 'none', 'MarkerFaceColor', 0.02);

                apply_generic(ax(i), 'xlim', [0.5 3.5], 'ylim', [-0.02 1.02], 'xticks', base,  ...
                    'yticks', [0 1],  ...
                    'tick_length', 0.03)
                xticklabels({'Cont', 'Shrt', 'Long'});
                axis square
                drawnow
            end
        end
    end

end


paths.fig = fullfile(paths.fig);

% if export
%     if quant.sm == 1
%         figure_title = sprintf('cmpr_pmix');
%         savefig(fh, figure_title)
%     else
%         figure_title = sprintf('cmpr_pmix');
%     end
%     
%     exporter(fh, paths, [figure_title, '.png'])
%     exporter(fh, paths, [figure_title, '.pdf'])
% end

