function [fh, ax, ax_inset] = fd_conditions(varargin)

opt = inputParser;
addParameter(opt, 'extra', []);
addParameter(opt, 'freezes', []);
addParameter(opt, 'results', []);
addParameter(opt, 'bin_size', 3);
addParameter(opt, 'conditions', false);
addParameter(opt, 'type', 'continuous');
addParameter(opt, 'no_y', true);
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
no_y = opt.Results.no_y;

col = cmapper();

dx = 1/60;
xbin = (0:dx:dx*(ceil(10/dx)+10))';
xxi     = 0:dx:1200;
h       = .05;

fs = 60;
bin_size_in_seconds = bin_size/fs;

if isempty(opt.Results.freezes)
    bouts = importdata(fullfile(opt.Results.results.bouts_path, 'surrogate.mat'));
    freezes = bouts;
end

results.fitted_model = sprintf('model_%s', results.fitted_model);
if ~isempty(results.points.censoring)
    freezes.durations_s(freezes.durations_s > results.points.censoring) = results.points.censoring + 1/60;
    freezes = freezes(freezes.durations_s >= results.points.truncation, :);

end

fh = figure('Position', [100 100 1250 580], 'Color', 'w');
t = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
% ax_inset = struct;
% ax = struct();
i = 0;

for idx_sm = 1:3
    for idx_ls = 1:2
        for idx_fs = 1:2

            i = i + 1;
            nexttile(t)
            hold on

            ax(i) = gca;

            [freezes_quant, ~] = quantilizer(freezes, 'idx_quanti', struct('sm', idx_sm, 'ls', idx_ls, 'fs', idx_fs));

            RTs{1,1} = freezes_quant.durations_s;
            RTs{2,1} = freezes_quant.durations_s;
            % RTD = kreg_single(RTs, RTs, xxi, xbin, h, min(freezes.durations_s) - 0.001, 500);

            censored = freezes_quant.durations_s > 10.5;
            [ks_noc] = ksdensity(freezes_quant.durations_s, xxi, 'BoundaryCorrection', 'reflection', 'Bandwidth', h, 'Support', [min(freezes.durations_s) - 0.001, max(freezes.durations_s) + 0.001 ]);
            [ks] = ksdensity(freezes_quant.durations_s, xxi, 'BoundaryCorrection', 'reflection', 'Bandwidth', h, 'Support', [min(freezes.durations_s) - 0.001, max(freezes.durations_s) + 0.001 ], 'Censoring', censored);

            censored_density = abs(trapz(xxi, ks) - trapz(xxi, ks_noc));

            % fill_between(xxi, RTD{1,2} + RTD{2,2}, RTD{1,2} - RTD{3,2}, [], 'FaceColor', 'r', 'FaceAlpha', 0.4, 'LineStyle', 'none')
            % plot(xxi, RTD{1,2}, 'Color', 'k', 'LineWidth', 2)

            histogram(freezes_quant.durations_s, min(freezes.durations_s) - 1/120:bin_size_in_seconds:600, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceColor',  col.empirical.(color)( 2*(idx_sm - 1) + idx_fs, :), 'FaceAlpha', 0.3)
            plot(xxi, ks, 'Color', col.empirical.(color)( 2*(idx_sm - 1) + idx_fs, :), 'LineWidth' , 2)

            %             hp = plot([-.3 -.3], [0 1], 'k-', 'LineWidth', 2, 'Clipping', 'off');
            %             text(-0.52, 0, '0', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 18)
            %             text(-0.52, 1, '1', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 18)

            
            if censored_inset

                inset_pos = ax(i).Position;  % Get the position of the current subplot
                inset_pos = [inset_pos(1) + 0.8 * inset_pos(3), inset_pos(2) + 0.56 * inset_pos(4), 0.015, 0.075];


                ax_inset(i) = axes('Position', inset_pos);  % Create inset axes with adjusted position
                hold(ax_inset(i), 'on');
                set(gca,'box','on')
                censored_x = results.points.censoring;
                scatter(censored_x, censored_density, 32, 'filled', 'MarkerFaceColor', col.empirical.(color)( 2*(idx_sm - 1) + idx_fs, :), 'MarkerEdgeColor', 'none');

                apply_generic(ax_inset(i), 'ylim', [0 0.3], 'xlim', [censored_x - 0.025 censored_x + 0.025], 'yticks', [0 0.3], 'xticks', censored_x, ...
                    'tick_length', 0.05, 'line_width', 2, 'no_x', true, 'font_size', 18);


            end

            axes(ax(i));
            scatter(censored_x, censored_density, 32, 'filled', 'MarkerFaceColor', col.empirical.(color)( 2*(idx_sm - 1) + idx_fs, :), 'MarkerEdgeColor', 'none');
            apply_generic(ax(i), 'xlim', [-0.1 10.6], 'ylim', [-0.05 2.05], 'no_y', no_y, 'font_size', 20, 'line_width', 2.2, 'yticks', [0 1], 'xticks', [results.points.truncation 10])
            drawnow
            
        end
    end
end

linkaxes(ax, 'xy');

if export
    paths.fig = results.fig_path;
    if conditions
        exporter(fh, paths, 'fits_xcondition.pdf')

    else
        exporter(fh, paths, 'fits.pdf')
    end
end
