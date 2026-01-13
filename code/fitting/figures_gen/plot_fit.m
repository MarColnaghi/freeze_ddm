function [fh, fd, f] = plot_fit(varargin)

opt = inputParser;
addParameter(opt, 'extra', []);
addParameter(opt, 'freezes', []);
addParameter(opt, 'results', []);
addParameter(opt, 'bin_size', 3);
addParameter(opt, 'conditions', false);
addParameter(opt, 'type', 'continuous');
addParameter(opt, 'censored_inset', false);
addParameter(opt, 'gt', false);
addParameter(opt, 'export', false);
addParameter(opt, 'paths', []);
addParameter(opt, 'no_dep', false);

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
no_dep = opt.Results.no_dep;

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
    if no_dep
        [~, f, fd] = nll_fly_ddm_newer(est_params, freezes(1,:), results.points, results.fitted_model, 'iid', 'p', extra);
    else
        [~, f, fd] = nll_fly_ddm_newer(est_params, freezes, results.points, results.fitted_model, 'iid', 'p', extra);

    end

    fh = figure('Position', [100 100 800 500], 'Color', 'w');
    hold on
    hh = histogram(freezes.durations_s, -1/120:bin_size_in_seconds:12, 'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none');
    
    % Robustly downsample fd
    fd_ds = downsample_vec(fd, bin_size, 'mean');

    if strcmp('discrete', type)
        f_ds = downsample_vec(f, bin_size, 'sum');
        f_ds = f_ds ./ bin_size_in_seconds;
        plot(fd_ds, f_ds, 'k--', 'LineWidth', 2)
        trapz(fd_ds(results.points.truncation <=  fd_ds), f_ds(results.points.truncation <= fd_ds))
    else
        f_ds = downsample_vec(f, bin_size, 'mean');
        % Handle the specific scaling for the final element if needed
        f_ds(end) = f_ds(end) ./ bin_size_in_seconds; 
        plot(fd_ds, f_ds, 'k--', 'LineWidth', 2)
        plot(fd, f, 'g--', 'LineWidth', 2)

    end
    
    if gt_plot 
        gt = table2array(results.ground_truth(1,~ismissing(results.estimates_mean)));
        [~, f_gt, fd_gt] = nll_fly_ddm_newer(gt, freezes, results.points, results.fitted_model, 'iid', 'p', extra);
        fd_ds_gt = downsample_vec(fd_gt, bin_size, 'mean'); 
        f_ds_gt = downsample_vec(f_gt, bin_size, 'sum');
        plot(fd_ds_gt, f_ds_gt ./ bin_size_in_seconds, 'b--', 'LineWidth', 2)
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
        xlim([censored_x - 0.1 censored_x + 0.1])
        scatter(results.points.censoring, f_ds(end), 240, '_', 'k', 'LineWidth', 2)
        xticks(results.points.censoring); xticklabels('cens'); xtickangle(0);
        apply_generic(ax_inset)
    end
else
    % ... [Code for 'conditions' = true] ...
    % Apply the same downsample_vec pattern inside your loop
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
                if isfield(extra, 'soc_mot_array')
                    ec.soc_mot_array = extra.soc_mot_array(mask, :);
                end
                [~, f, fd] = nll_fly_ddm_newer(est_params, freezes_quant, results.points, results.fitted_model, 'iid', 'p', ec);
                histogram(freezes_quant.durations_s, 1/120:bin_size/fs:12, 'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none')
                
                fd_ds = downsample_vec(fd, bin_size, 'mean');
                if strcmp('discrete', type)
                    f_ds = downsample_vec(f, bin_size, 'sum');
                    f_ds = f_ds ./ bin_size_in_seconds;
                else
                    f_ds = downsample_vec(f, bin_size, 'mean');
                    f_ds(end) = f_ds(end) ./ bin_size_in_seconds;
                end
                plot(fd_ds, f_ds, 'k--', 'LineWidth', 2)
                
                % [Inset logic remains the same using updated f_ds]
                ax(i) = gca;
                apply_generic(gca); xlim([-0.1 11.1]); ylim([-0.001 0.801]); yticks([]);
                inset_pos = ax(i).Position;
                inset_pos = [inset_pos(1) + 0.6 * inset_pos(3), inset_pos(2) + 0.4 * inset_pos(4), 0.05, 0.075];
                ax_inset(i) = axes('Position', inset_pos); hold on;
                set(gca,'box','on')
                h_ins_h = histogram(freezes_quant.durations_s, results.points.censoring + 0.002:bin_size_in_seconds:12, 'Normalization', 'pdf', 'FaceColor', 'r', 'EdgeColor', 'none');
                [~, idx] = find(h_ins_h.Values > 0);
                if any(idx), height(i) = h_ins_h.Values(idx); else height(i) = nan; end
                xlim([results.points.censoring - 0.1 results.points.censoring + 0.1])
                scatter(results.points.censoring, f_ds(end), 240, '_', 'k', 'LineWidth', 2)
                apply_generic(ax_inset(i)); hold off;
            end
        end
    end
    linkaxes(ax); linkaxes(ax_inset);
    try ylim(ax_inset, [0; max(height) + 1.5]); catch; end
end

if export
    paths.fig = results.fig_path;
    exporter(fh, paths, iif(conditions, 'fits_xcondition.pdf', 'fits.pdf'));
end
end

% Helper Function for Shape-Independent Downsampling
function out = downsample_vec(v, bin_size, mode)
    v = v(:)'; % Always work with a row vector
    n = numel(v);
    num_bins = floor((n - 1) / bin_size);
    main_part = v(1:num_bins * bin_size);
    reshaped_part = reshape(main_part, bin_size, []);
    
    if strcmpi(mode, 'sum')
        reduced = sum(reshaped_part, 1);
    else
        reduced = mean(reshaped_part, 1);
    end
    
    out = [reduced, v(end)]; % Append the last element (censored bin)
end