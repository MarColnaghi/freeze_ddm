function [fh, fd, f] = plot_fit(varargin)

opt = inputParser;
addParameter(opt, 'extra', []);
addParameter(opt, 'freezes', []);
addParameter(opt, 'results', []);
addParameter(opt, 'bin_size', 5);
addParameter(opt, 'conditions', false);
addParameter(opt, 'type', 'continuous');

addParameter(opt, 'censored_inset', false);
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

fs = 60;
bin_size_in_seconds = bin_size/fs;

if isempty(opt.Results.freezes)
    load(fullfile(opt.Results.results.bouts_path, 'surrogate.mat'));
    freezes = surrogate;
end
if isempty(opt.Results.extra)
    load(fullfile(opt.Results.results.bouts_path, 'extra.mat'));
end

est_params = table2array(results.estimates_mean(:, ~ismissing(results.estimates_mean)));

results.fitted_model = sprintf('model_%s', results.fitted_model);
if ~isempty(results.points.censoring)
    freezes.durations_s(freezes.durations_s > results.points.censoring) = results.points.censoring + 1/60;
end

if ~conditions

    [~, f, fd] = nll_fly_ddm_newer(est_params, freezes, results.points, results.fitted_model, 'iid', 'p', extra);
    fh = figure('Position', [100 100 800 500], 'Color', 'w');
    hold on
    histogram(freezes.durations_s, 1/120:bin_size_in_seconds:12, 'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none')
    fd_ds = [mean(reshape(fd(1:end - 1), bin_size, []), 1) fd(end)]; % f_ds = [sum(reshape(f(1:end - 1), bin_size, []), 1) f(end)];
    %plot(fd_ds, f_ds , 'k--', 'LineWidth', 2)

    if strcmp('discrete', type)
        f_ds = [sum(reshape(f(1:end - 1), bin_size, []), 1) f(end)];
        f_ds = f_ds ./ bin_size_in_seconds;
        plot(fd_ds, f_ds, 'k--', 'LineWidth', 2)

    else
        f_ds = [mean(reshape(f(1:end - 1), bin_size, []), 1) f(end) ./ bin_size_in_seconds];
        plot(fd_ds, f_ds, 'k--', 'LineWidth', 2)

    end

    if gt_plot % This has to be fixed
        gt = table2array(results.ground_truth(1,~ismissing(results.estimates_mean)));

        [~, f, fd] = nll_fly_ddm_newer(gt, freezes, results.points, results.fitted_model, 'iid', 'p', extra);
        fd_ds = [mean(reshape(fd(1:end - 1), bin_size, []), 1) fd(end)]; f_ds = [sum(reshape(f(1:end - 1), bin_size, []), 1) f(end)];
        plot(fd_ds, f_ds ./ bin_size_in_seconds, 'b--', 'LineWidth', 2)
    end

    apply_generic(gca)
    xlabel('Freeze Duration (s)')
    ylabel('pdf')
    xlim([-0.1 11.1])
    ylim([-0.001 0.501])

    if censored_inset
        ax_inset = axes('Position', [0.6 0.5 0.1 0.3]); hold on;
        h_ins_h = histogram(freezes.durations_s, results.points.censoring+1e-2:bin_size_in_seconds:12, 'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none');
        [~, idx] = find(h_ins_h.Values > 0);
        censored_x = results.points.censoring;

        if any(idx)
            height = h_ins_h.Values(idx);
        else
            height = nan;
        end
        xlim([censored_x - 0.1 censored_x + 0.1])
        scatter(results.points.censoring, f_ds(end), 240, '_', 'k', 'LineWidth', 2)
        xticks(results.points.censoring); xticklabels('cens'); xtickangle(0);
        apply_generic(ax_inset)
    end

else
    fh = figure('Position', [100 100 1400 900], 'Color', 'w');
    t = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'loose');
    i = 0;
    for idx_sm = 1:3
        for idx_ls = 1:2
            for idx_fs = 1:2

                i = i + 1;
                nexttile(t)
                hold on

                [freezes_quant, mask] = quantilizer(freezes, 'idx_quanti', struct('sm', idx_sm, 'ls', idx_ls, 'fs', idx_fs));
                ec = extra;
                ec.soc_mot_array = extra.soc_mot_array(mask, :);
                [~, f, fd] = nll_fly_ddm_newer(est_params, freezes_quant, results.points, results.fitted_model, 'iid', 'p', ec);

                histogram(freezes_quant.durations_s, 1/120:bin_size/fs:12, 'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none')
                fd_ds = [mean(reshape(fd(1:end - 1), bin_size, []), 1) fd(end)];

                if strcmp('discrete', type)
                    f_ds = [sum(reshape(f(1:end - 1), bin_size, []), 1) f(end)];
                    f_ds = f_ds ./ bin_size_in_seconds;
                    plot(fd_ds, f_ds, 'k--', 'LineWidth', 2)

                else
                    f_ds = [mean(reshape(f(1:end - 1), bin_size, []), 1) f(end) ./ bin_size_in_seconds];
                    plot(fd_ds, f_ds, 'k--', 'LineWidth', 2)

                end

                %
                %                 est_params_first = est_params;
                %                 est_params_first([7 8 9 10]) = 10 * ones(1,4);
                %                 [~, f, fd] = nll_fly_ddm_newer(est_params_first, freezes_quant, results.points, results.fitted_model, 'iid', 'p', ec);
                %                 histogram(freezes_quant.durations_s, results.points.truncation:bin_size/fs:12, 'Normalization', 'probability', 'FaceColor', 'r', 'EdgeColor', 'none')
                %                 fd_ds = [mean(reshape(fd(1:end - 1), bin_size, []), 1) fd(end)]; f_ds = [sum(reshape(f(1:end -1), bin_size, []), 1) f(end)];
                %                 plot(fd_ds, f_ds, 'g--', 'LineWidth', 1)
                %
                %                 est_params_second = est_params;
                %                 est_params_second([7 8 9 10]) = - 10 * ones(1,4);
                %                 [~, f, fd] = nll_fly_ddm_newer(est_params_second, freezes_quant, results.points, results.fitted_model, 'iid', 'p', ec);
                %                 histogram(freezes_quant.durations_s, results.points.truncation:bin_size/fs:12, 'Normalization', 'probability', 'FaceColor', 'r', 'EdgeColor', 'none')
                %                 fd_ds = [mean(reshape(fd(1:end - 1), bin_size, []), 1) fd(end)]; f_ds = [sum(reshape(f(1:end -1), bin_size, []), 1) f(end)];
                %                 plot(fd_ds, f_ds, 'b--', 'LineWidth', 1)


                
                ax(i) = gca;
                apply_generic(gca)
                xlim([-0.1 11.1])
                ylim([-0.001 0.801])
                yticks([])

                inset_pos = ax(i).Position;  % Get the position of the current subplot
                inset_pos = [inset_pos(1) + 0.6 * inset_pos(3), inset_pos(2) + 0.4 * inset_pos(4), 0.05, 0.075];
                ax_inset(i) = axes('Position', inset_pos);  % Create inset axes with adjusted position
                hold(ax_inset(i), 'on');
                set(gca,'box','on')
                h_ins_h = histogram(freezes_quant.durations_s, results.points.censoring + 0.002:bin_size_in_seconds:12, 'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none');
                [~, idx] = find(h_ins_h.Values > 0);
                censored_x = results.points.censoring;

                if any(idx)
                    height(i) = h_ins_h.Values(idx);
                else
                    height(i) = nan;
                end
                xlim([censored_x - 0.1 censored_x + 0.1])
                scatter(results.points.censoring, f_ds(end), 240, '_', 'k', 'LineWidth', 2)
                xticks(results.points.censoring); xticklabels('cens'); xtickangle(0);
                apply_generic(ax_inset(i))
                hold(ax_inset(i), 'off');
                

            end
            drawnow
        end
    end
    linkaxes(ax)
    linkaxes(ax_inset)
    ylim([0; max(height) + 1.5])
end
if export
    paths.fig = results.fig_path;
    if conditions
        exporter(fh, paths, 'fits_xcondition.pdf')

    else
        exporter(fh, paths, 'fits.pdf')
    end
end
