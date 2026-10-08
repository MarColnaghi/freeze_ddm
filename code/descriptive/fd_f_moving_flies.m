id_code = 'imm3_mob3_pc4';
paths_out = path_generator('folder', 'experimental_vars/freeze_durations', 'bouts_id', id_code, 'imfirst', false);

col = cmapper();
thresholds = define_thresholds;
thresholds.le_window_fl = [5 40];
thresholds.le_window_sl = [15 50];

bouts = importdata(fullfile(paths_out.dataset, 'bouts.mat'));
bouts = bouts_formatting(bouts, thresholds);

trunc_point = 60;
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', trunc_point);
threshold = 70;
[bouts_proc, contact_mask] = impose_contact_threshold(bouts_proc, 'threshold', threshold);

bouts_proc.ends = bouts_proc.onsets + bouts_proc.ending_time;
bouts_proc.durations_s = bouts_proc.ending_time ./ 60;
bouts_proc.durations = bouts_proc.ending_time;

sm_freeze_full = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat','norm_factor', 10, 'cache', 'motion_cache');
window_dur = 60;
sm_freeze_ili = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat','norm_factor', 10, 'cache', 'motion_cache', 'window', [0 window_dur]);

bouts_proc.sm = mean(sm_freeze_full, 2, 'omitnan');
%bouts_proc.sm = mean(sm_freeze_full, 2, 'omitnan');

bouts_proc = bouts_proc(bouts_proc.durations > trunc_point, :);
bouts_proc = bouts_proc(bouts_proc.sm < 1.5,:);

fh = figure('color', 'w', 'Position', [100, 100, 350, 600]);

tl = tiledlayout(5, 1, 'TileSpacing', 'none', 'Padding', 'none');


%% Freezing-bout duration densities by number of moving flies
xmax   = 10;                              % right edge of plotted density (s)
gap    = 0.2;                             % break between curve and bar (s)
h      = 0.12;                            % kernel bandwidth (s)
ymax   = 1;                               % shared y-limit
lo     = min(bouts_proc.durations_s);     % reflection point (detection floor)
groups = 4:-1:0;
n      = numel(groups);

% ---- pass 1: KM weights, density, tail mass ----------------------------
S  = zeros(1, n);
xg = cell(1, n);
fg = cell(1, n);
fh = figure('color', 'w', 'Position', [100, 100, 400, 600]);

for k = 1:n
    b = bouts_proc(bouts_proc.moving_flies == groups(k), :);

    if isempty(b) || all(b.is_censored)   % no uncensored events to smooth
        S(k) = 1;  xg{k} = [];  fg{k} = [];
        continue
    end

    [F, t] = ecdf(b.durations_s, 'Censoring', b.is_censored);

    i10  = find(t <= xmax, 1, 'last');
    S(k) = 1 - F(i10);                    % KM survivor at xmax

    w  = diff(F);                         % KM jump sizes (column)
    te = t(2:end).';                      % uncensored event times (row)

    xk    = (lo : 1/60 : min(xmax, te(end))).';
    fg{k} = (normpdf((xk - te)/h) + normpdf((xk + te - 2*lo)/h)) * w / h;
    xg{k} = xk;
end

% bar width chosen so the tallest bar clears the y-limit
bw = max(0.6, max(S) / (0.9 * ymax));
xb = xmax + gap;

% ---- pass 2: plot ------------------------------------------------------
tl = tiledlayout(n, 1, 'TileSpacing', 'none', 'Padding', 'none');

for k = 1:n
    nexttile
    hold on
    c = col.vars.moving_flies(groups(k) + 1, :);
    x = xg{k};
    f = fg{k};

    if ~isempty(x)
        fill([x; flipud(x)], [f; zeros(size(f))], c, ...
             'EdgeColor', 'none', 'FaceAlpha', 0.3)
        plot(x, f, 'Color', c, 'LineWidth', 2.5)
    end

    fill([xb xb+bw xb+bw xb], [0 0 S(k)/bw S(k)/bw], c, ...
         'FaceAlpha', 0.3, 'EdgeColor', c, 'LineWidth', 1.5)

    if groups(k) == 0
        apply_generic(gca, 'ylim', [-0.05 ymax], 'xlim', [-0.1, 11], ...
                      'no_y', true, 'tick_length', 0.02, 'font_size', 30)
        xticks([0 5 10])

    else
        apply_generic(gca, 'ylim', [-0.05 ymax], 'xlim', [-0.1, 11], ...
                      'no_y', true, 'tick_length', 0.02, 'no_x', true, ...
                      'font_size', 30)
    end

    
    if ~isempty(x)
        fprintf('%d moving: n=%4d, censored=%3.0f%%, area=%.3f, S(>%g)=%.3f, sum=%.3f\n', ...
                groups(k), height(bouts_proc(bouts_proc.moving_flies == groups(k), :)), ...
                100*mean(bouts_proc.is_censored(bouts_proc.moving_flies == groups(k))), ...
                trapz(x, f), xmax, S(k), trapz(x, f) + S(k));
    end
end
exporter(fh, paths, 'fd_f_moving_flies.pdf')
%%
fh = figure('color', 'w', 'Position', [100, 100, 500, 500]);

tl = tiledlayout(1, 1, 'TileSpacing', 'none', 'Padding', 'none');


for idx_moving_flies = 4:-1:0
    hold on
    bouts_quant = bouts_proc(bouts_proc.moving_flies == idx_moving_flies, :);
    [fkde, x] = ksdensity(bouts_quant.durations_s, 0:1/60:max(bouts_proc.durations_s), ...
        'BoundaryCorrection', 'reflection', 'Bandwidth', 0.12, 'Censoring', bouts_quant.is_censored,...
        'Support', [min(bouts_proc.durations_s) - 1e-12, max(bouts_proc.durations_s) + 1e-12]);
    plot(x, fkde, 'Color', col.vars.moving_flies(idx_moving_flies + 1, :), 'LineWidth', 1.8)
    if idx_moving_flies == 0
        apply_generic(gca, 'ylim', [0 1.3], 'xlim', [0 10], 'no_y', true, 'tick_length', 0.02, 'no_x', false)
    else
        apply_generic(gca, 'ylim', [0 1.3], 'xlim', [0 10], 'no_y', true, 'tick_length', 0.02, 'no_x', true)
    end


end


%%

%%
fh = figure('color', 'w', 'Position', [100, 100, 500, 500]);

tl = tiledlayout(1, 1, 'TileSpacing', 'none', 'Padding', 'none');

for idx_sm = 1:4
    hold on
    bouts_quant = quantilizer_v2(bouts_proc, 'total_quantiles', struct('sm', 4, 'fs', 1), 'indexed_quantile', struct('sm', idx_sm, 'fs', 1));
    [fkde, x] = ksdensity(bouts_quant.durations_s, 0:1/60:max(bouts_proc.durations_s),'Censoring', bouts_quant.is_censored, ...
        'BoundaryCorrection', 'reflection', 'Bandwidth', 0.12, ... 
        'Support', [min(bouts_proc.durations_s) - 1e-12, max(bouts_proc.durations_s) + 1e-12], 'Function', 'pdf');
    plot(x, fkde, 'Color', col.vars.sm(idx_sm, :), 'LineWidth', 1.8)
        apply_generic(gca, 'ylim', [0 1.3], 'xlim', [0 10], 'no_y', true, 'tick_length', 0.02, 'no_x', false)



end