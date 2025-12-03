clearvars
paths = path_generator('folder', '/social_motion');
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));
pixel_cache = importdata(fullfile(paths.cache_path, 'pixel_cache.mat'));
loom_cache = importdata(fullfile(paths.cache_path, 'loom_cache.mat'));
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
thresholds = define_thresholds;

window = [-601 630];
offsets = (window(1) : window(2));
total_looms = 20; 
total_flies = size(motion_cache, 1);

n_moving_flies = nan(total_flies, 1);
sloom = nan(total_flies, 1);
sm_around_loom = nan(total_flies, total_looms, length(offsets));

for idx_fly = 1:total_flies
    sm_fly = motion_cache(idx_fly);
    loom_frames = find(diff(loom_cache(idx_fly)) == 1);
    idx_slice = loom_frames(:) + offsets;
    fly_loom_x_sm = sm_fly(idx_slice);
    n_moving_flies(idx_fly) = unique(bouts.moving_flies(bouts.fly == idx_fly));
    sloom(idx_fly) = unique(bouts.sloom(bouts.fly == idx_fly));

    sm_around_loom(idx_fly, :, :) = fly_loom_x_sm;

end

col = cmapper( '', 30);

%% Plot 4 scatters for the different windows

bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'min_dur', 30, 'sloom', 25);
bouts_temp = bouts_proc;

avg_sm_window = table();
windows = [-600 -300 -180 -60];
window = -1200;
avg_sm_freeze = nan(height(bouts_temp), 1);
avg_sm_600 = nan(height(bouts_temp), 1);
avg_sm_300 = nan(height(bouts_temp), 1);
avg_sm_180 = nan(height(bouts_temp), 1);
avg_sm_60 = nan(height(bouts_temp), 1);

for idx_bouts = 1:height(bouts_temp)
    fly_id = bouts_temp.fly(idx_bouts);
    sm_fly = motion_cache(fly_id);

    ons = bouts_temp.onsets(idx_bouts);
    off = bouts_temp.ends(idx_bouts) - 1;
    sm_freeze = sm_fly(ons:off);
    avg_sm_freeze(idx_bouts) = mean(sm_freeze);
    
    sm_b4_freeze = sm_fly(ons + window(1):ons - 1);

    avg_sm_600(idx_bouts) = mean(sm_b4_freeze(-window(1) + windows(1): (- window(1))));
    avg_sm_300(idx_bouts) = mean(sm_b4_freeze(-window(1) + windows(2): (- window(1))));
    avg_sm_180(idx_bouts) = mean(sm_b4_freeze(-window(1) + windows(3): (- window(1))));
    avg_sm_60(idx_bouts) =  mean(sm_b4_freeze(-window(1) + windows(4): (- window(1))));

end

avg_sm_window.avg_sm_60 = avg_sm_60;
avg_sm_window.avg_sm_180 = avg_sm_180;
avg_sm_window.avg_sm_300 = avg_sm_300;
avg_sm_window.avg_sm_600 = avg_sm_600;
avg_sm_window.avg_sm_freeze = avg_sm_freeze;
 
col = cmapper( '', 30);
fh = figure('color','w','Position',[100, 100, 1300, 500]);
tiledlayout(1, 4, 'TileSpacing', 'loose', 'Padding', 'loose')

deltas = nan(4, height(bouts_temp));
size = 12;

for idx_window = 1:4
    
    nt = nexttile;

    hold on
    plot([-10 30], [-10 30], 'k--', 'LineWidth', 1)
    current_window = windows(idx_window);
    scatter(bouts_temp.avg_sm, avg_sm_window.(sprintf('avg_sm_%d', -current_window)),  size, bouts_temp.durations_s, 'o', 'filled', 'MarkerFaceAlpha', .35)
    deltas(idx_window, :) = avg_sm_window.(sprintf('avg_sm_%d', -current_window)) - bouts_temp.avg_sm;
    axis square
    apply_generic(gca, 'xlim', [-1 21], 'ylim', [-1 21], 'xticks', [0 20], 'yticks', [0 20])
    clim([0 4])

    xlabel({'Avg. Soc.', 'Mot. Freeze'}, 'FontSize', 20)
    ylabel({'Avg. Soc.', 'Mot. Before'}, 'FontSize', 20)

    ax = gca;
    txtpos = ax.Position;
    txtpos(2) = txtpos(2) + 0.055;
    annotation('textbox', txtpos, 'String', sprintf('Window: [%d 0]', current_window) ,'HorizontalAlignment', 'left', 'VerticalAlignment','top', ...
        'LineStyle', 'none', 'FontSize', 18);
end

colormap(colorcet('L20'))
exporter(fh, paths, 'scatters_sm_pre_post.pdf')

%% Plot the corrcoeff


fh = figure('color','w','Position', [100, 100, 600, 350]);
tiledlayout(1, 1, 'TileSpacing', 'compact')
nexttile
ax = gca;
hold on

for idx_sl = [unique(bouts.sloom)]'

    bouts_proc = data_parser_new(bouts, 'type', 'immobility', ...
        'period', 'loom', 'window', 'le', 'min_dur', 360, 'nloom', 1:20);
    bouts_proc = bouts_proc(bouts_proc.sloom == idx_sl, :);

  %  for idx_mov_flies = 0:4
        %idx_mov_flies = 0;
   bouts_temp = bouts_proc(bouts_proc.moving_flies == 0, :);

        max_window = 1200;
        n_bouts = height(bouts_temp);

        avg_sm_freeze = nan(n_bouts, 1);
        avg_sm_pre    = nan(n_bouts, max_window);

        for idx_bouts = 1:n_bouts
            fly_id = bouts_temp.fly(idx_bouts);
            sm_fly = motion_cache(fly_id);

            ons = bouts_temp.onsets(idx_bouts);
            off = bouts_temp.ends(idx_bouts) - 1;

            sm_freeze = sm_fly(ons:off);
            avg_sm_freeze(idx_bouts) = mean(sm_freeze);

            start_idx = max(1, ons - max_window);
            sm_b4_freeze = sm_fly(start_idx:ons - 1);

            L = numel(sm_b4_freeze);
            sm_rev = sm_b4_freeze(end:-1:1);

            cs = cumsum(sm_rev);
            means_all = cs ./ (1:L)';
            avg_sm_pre(idx_bouts, 1:L) = means_all;
        end


        coeffs = corrcoef([avg_sm_freeze avg_sm_pre]);
        plot(-max_window:1:-1, coeffs(1, end:-1:2), 'LineWidth', 2, 'Color', col.loomspeed(1 + idx_sl / 25, :))

   % end
end

xline(-60, 'k--')
apply_generic(ax, 'ylim', [0.0 0.8], 'xlim', [-600 0], 'xticks', sort([-600 0]))

%%

fh = figure('color','w','Position', [100, 100, 600, 350]);
tiledlayout(1, 1, 'TileSpacing', 'compact')
nexttile
ax = gca;
hold on

for idx_sl = [unique(bouts_proc.sloom_norm)]'

    bouts_proc = data_parser_new(bouts, 'type', 'immobility', ...
        'period', 'loom', 'window', 'le', 'min_dur', 30, 'nloom', 1:20);
    bouts_temp = bouts_proc(bouts_proc.sloom_norm == idx_sl, :);

    max_window = 1200;
    n_bouts = height(bouts_temp);

    avg_sm_freeze = nan(n_bouts, 1);
    avg_sm_pre    = nan(n_bouts, max_window);

    for idx_bouts = 1:n_bouts
        fly_id = bouts_temp.fly(idx_bouts);
        sm_fly = motion_cache(fly_id);

        ons = bouts_temp.onsets(idx_bouts);
        off = bouts_temp.ends(idx_bouts) - 1;
        loom_ons = bouts_temp.loom_ts(idx_bouts);

        sm_freeze = sm_fly(ons:off);
        avg_sm_freeze(idx_bouts) = mean(sm_freeze);

        start_idx = max(1, loom_ons - max_window);
        sm_b4_freeze = sm_fly(start_idx:loom_ons - 1);

        L = numel(sm_b4_freeze);
        sm_rev = sm_b4_freeze(end:-1:1);

        cs = cumsum(sm_rev);
        means_all = cs ./ (1:L)';
        avg_sm_pre(idx_bouts, 1:L) = means_all;
    end


    coeffs = corrcoef([avg_sm_freeze avg_sm_pre]);
    plot(-max_window:1:-1, coeffs(1, end:-1:2), 'LineWidth', 2, 'Color', col.loomspeed(1 + idx_sl, :))

end

xline(-60, 'k--')
apply_generic(ax, 'ylim', [0.3 0.7], 'xlim', [-600 0], 'xticks', sort([-600 0]))



%%
fh = figure('color','w','Position',[100, 100, 900, 400]);
tiledlayout(1, 5, 'TileSpacing', 'compact', 'Padding', 'loose')

for idx_mov_flies = 0:4
    nexttile
    mask = bouts_temp.moving_flies == idx_mov_flies;
    histogram(deltas(4, mask), -30:.5:30, 'Normalization', 'count', 'Orientation', 'horizontal')
    hold on
    apply_generic(gca, 'xlim', [0 100], 'ylim', [-30 30])
end


%% 
fh = figure('color','w','Position',[100, 100, 900, 400]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'loose')
nexttile
scatter(deltas(4,:), bouts_temp.onsets_loomaligned)
nexttile
scatter(deltas(4,:), bouts_temp.durations_s)
nexttile
scatter(deltas(4,:), bouts_temp.nloom)


%%
[sorted_deltas, idx] = sort(deltas(1,:));
figure
hold on
for idx_bouts = [idx(1)]
    fly_id = bouts_temp.fly(idx_bouts);
    sm_fly = motion_cache(fly_id);

    ons = bouts_temp.onsets(idx_bouts) - 300;
    off = bouts_temp.ends(idx_bouts) - 1;
    sm_freeze = sm_fly(ons:off);
    plot(sm_freeze)
end