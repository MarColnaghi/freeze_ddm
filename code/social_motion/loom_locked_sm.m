tic
clearvars
paths = path_generator('folder', '/social_motion');
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));
pixel_cache = importdata(fullfile(paths.cache_path, 'pixel_cache.mat'));
loom_cache = importdata(fullfile(paths.cache_path, 'loom_cache.mat'));
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
thresholds = define_thresholds;

window = [-300 300];
offsets = (window(1) : window(2));
total_looms = 20;
total_flies = size(motion_cache, 1);

n_moving_flies = nan(total_flies, 1);
sloom = nan(total_flies, 1);
sm_around_loom = nan(total_flies, total_looms, length(offsets));

for idx_fly = 1:total_flies
    sm_fly = motion_cache(idx_fly);
    freeze_frames = find(diff(loom_cache(idx_fly)) == 1);
    idx_slice = freeze_frames(:) + offsets;
    fly_loom_x_sm = sm_fly(idx_slice);
    n_moving_flies(idx_fly) = unique(bouts.moving_flies(bouts.fly == idx_fly));
    sloom(idx_fly) = unique(bouts.sloom(bouts.fly == idx_fly));

    sm_around_loom(idx_fly, :, :) = fly_loom_x_sm;

end

% Plot Time Locked Social Motion

col = cmapper( '', 30);
fh = figure('color','w','Position',[100, 100, 1400, 630]);
tiledlayout(2, 5, 'TileSpacing', 'compact')

graphvector.yaxis = [1 0 0 0 0 0 0 0 0 0];
y_pos = -2; y_pos_scale = 0.75;

for idx_sloom = unique(sloom)'

    for idx_moving_flies = unique(n_moving_flies)'

        nt = nexttile;
        hold on

        mask = n_moving_flies == idx_moving_flies & sloom == idx_sloom;
        sm_loom = squeeze(mean(sm_around_loom(mask, :, :), 1));

        colororder(col.vars.ln(end - total_looms:end, :))
        plot(offsets, sm_loom')



        colororder(col.vars.ln(end - total_looms:end, :))

        apply_generic(gca, 'no_y', true, 'ylims', [-1 21], 'xlims', [-60 180], ...
      'xticks', [-60 0 180], 'yticks', [0 10 20], 'tick_length', 0.04)

        xl = xlim;   yl = ylim;
        wx = diff(xl);
        wy = diff(yl);

        plot([xl(1) - 0.02*wx xl(1) - 0.02*wx], [0 10], 'k-', 'LineWidth', 2, 'Clipping', 'off');
       
        if ~(idx_moving_flies == 0 && idx_sloom == 25)
        else
            text(xl(1) - 0.04*wx, 0, '0', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Clipping', 'off', 'FontSize', 16)
            text(xl(1) - 0.04*wx, 10, '10', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Clipping', 'off', 'FontSize', 16)

        end

        % ---- BLACK SCALE BAR (time) ----
        x1 = 0;
        x2 = 30;

        y_top_k = yl(1) - 0.02*wy;
        y_bot_k = yl(1) - 0.05*wy;

        % ---- BLUE THRESHOLD BAR ----
        y_bot_b = yl(1) - 0.09*wy;
        y_top_b = yl(1) - 0.06*wy;

        if idx_sloom == 25
            x_thr = thresholds.le_window_sl;
        elseif idx_sloom == 50
            x_thr = thresholds.le_window_fl;
        else
            x_thr = [];
        end

        if idx_sloom == 25
            fill([x1 x1 x2 x2], ...
                [y_top_k y_bot_k y_bot_k y_top_k], ...
                'k', 'Clipping','off', 'EdgeColor','none');
            fill([x_thr fliplr(x_thr)], ...
                [y_bot_b y_bot_b y_top_b y_top_b], ...
                'b', 'Clipping','off', 'EdgeColor','none');
        elseif idx_sloom == 50
            fill([x1 x1 x2 x2], ...
                [y_top_k y_bot_k y_bot_k y_top_k], ...
                'k', 'Clipping','off', 'EdgeColor','none');
            fill([x_thr fliplr(x_thr)], ...
                [y_bot_b y_bot_b y_top_b y_top_b], ...
                'b', 'Clipping','off', 'EdgeColor','none');
        end

        if idx_sloom == 25
            xlabels = true;
            txtpos = nt.Position;
            txtpos(2) = txtpos(2) + 0.03;
            annotation('textbox', txtpos, 'String', sprintf('Moving Flies: %d', idx_moving_flies) ,'HorizontalAlignment', 'left', 'VerticalAlignment','top', ...
                'LineStyle', 'none', 'FontSize', 18);

        end

    end
end

exporter(fh, paths, 'loom_evoked_sm.pdf')

toc



%%


graphvector.yaxis = [1 0 0 0 0 0 0 0 0 0];

sloom = sloom(~isnan(sloom));
n_moving_flies = n_moving_flies(~isnan(n_moving_flies));

for idx_sloom = unique(sloom)'

    for idx_moving_flies = unique(n_moving_flies)'

        fh = figure('color','w','Position',[100, 100, 800, 1080]);
        tiledlayout(4, 5, 'TileSpacing', 'compact')

        mask = n_moving_flies == idx_moving_flies & sloom == idx_sloom;
        sm_loom_mean = squeeze(mean(sm_around_loom(mask, :, :), 1, 'omitnan'));
        sm_loom_std = squeeze(std(sm_around_loom(mask, :, :), 1, 'omitnan'));


        for idx_nloom = 1:20

            nt = nexttile;
            hold on

            curr_loom_mean = sm_loom_mean(idx_nloom, :);
            curr_loom_std = sm_loom_std(idx_nloom, :);

            plot(offsets, squeeze(sm_around_loom(mask, idx_nloom, :)), 'Color', [col.vars.ln(10 + idx_nloom, :) 0.1], 'LineWidth', 0.02)
            plot(offsets, curr_loom_mean, 'k-', 'LineWidth', 2)

            %[hfb] = fill_between(offsets, curr_loom_mean - curr_loom_std, curr_loom_mean + curr_loom_std, []);
            %hfb.EdgeColor = 'none'; hfb.FaceAlpha = 0.3;
            apply_generic(gca, 'no_y', true, 'ylims', [-2 22], 'xlims', [-60 180], ...
                'no_xlabels', xlabels, 'xticks', [-60 0 180], 'tick_length', 0.04)
        end


    end
end







%%
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'sloom', 50);
quantile_sm = quantile(bouts_proc.avg_sm_freeze_norm, 0/10);
bouts_temp = bouts_proc(bouts_proc.avg_sm_freeze_norm >= quantile_sm, :);
figure
hold on
for idx_bouts = 1:height(bouts_temp)
    fly_id = bouts_temp.fly(idx_bouts);
    sm_fly = motion_cache(fly_id);

    ons = bouts_temp.onsets(idx_bouts);
    off = bouts_temp.ends(idx_bouts) - 1;
    sm_freeze = sm_fly(ons:off);
    z(idx_bouts) = mean(sm_freeze);
    %plot(-length(sm_freeze) + 1:0, sm_freeze)

end

figure
histogram(bouts_temp.durations_s)

%% Plot the freezes

figure
bouts_proc = bouts(bouts.moving_flies == 3, :);
bouts_proc = data_parser_new(bouts_proc, 'type', 'immobility', 'period', 'loom', 'window', 'le');
bouts_proc = bouts_proc(bouts_proc.avg_sm_freeze_norm < 0.3, :);

scatter3(bouts_proc.durations, bouts_proc.nloom_norm, bouts_proc.avg_sm_freeze_norm, '.')
hold on
bouts_proc = bouts(bouts.moving_flies == 3, :);

bouts_proc = data_parser_new(bouts_proc, 'type', 'immobility', 'period', 'loom', 'window', 'le');
bouts_proc = bouts_proc(bouts_proc.avg_sm_freeze_norm > 0.3, :);

scatter3(bouts_proc.durations, bouts_proc.nloom_norm, bouts_proc.avg_sm_freeze_norm, '.')

xlim([0 1200])

%%

