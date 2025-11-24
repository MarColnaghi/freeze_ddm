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


%% Plot different correlations

col = cmapper( '', 30);
fh = figure('color','w','Position',[100, 100, 1400, 400]);

windows = [-600 -300 -180 -60];
tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'loose')


for idx_window = 1:4
    current_window = windows(idx_window);
    c = mean(sm_around_loom(:, :, -window(1) + current_window: (- window(1))), 3);
    scatter()

end

%%

bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'min_dur', 30, 'sloom', 25);
quantile_sm = quantile(bouts_proc.avg_sm_freeze_norm, 0/10);
bouts_temp = bouts_proc(bouts_proc.avg_sm_freeze_norm >= quantile_sm, :);

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
    
    sm_freeze = sm_fly(ons + window(1):ons - 1);

    avg_sm_600(idx_bouts) = mean(sm_freeze(-window(1) + windows(1): (- window(1))));
    avg_sm_300(idx_bouts) = mean(sm_freeze(-window(1) + windows(2): (- window(1))));
    avg_sm_180(idx_bouts) = mean(sm_freeze(-window(1) + windows(3): (- window(1))));
    avg_sm_60(idx_bouts) =  mean(sm_freeze(-window(1) + windows(4): (- window(1))));

    %plot(-length(sm_freeze) + 1:0, sm_freeze)

end
avg_sm_window.avg_sm_60 = avg_sm_60;
avg_sm_window.avg_sm_180 = avg_sm_180;
avg_sm_window.avg_sm_300 = avg_sm_300;
avg_sm_window.avg_sm_600 = avg_sm_600;
avg_sm_window.avg_sm_freeze = avg_sm_freeze;
 
col = cmapper( '', 30);
fh = figure('color','w','Position',[100, 100, 1400, 400]);

tiledlayout(1, 5, 'TileSpacing', 'compact', 'Padding', 'loose')

deltas = nan(4, height(bouts_temp));
for idx_window = 1:4
    nexttile
    current_window = windows(idx_window);
    scatter(bouts_temp.avg_sm, avg_sm_window.(sprintf('avg_sm_%d', -current_window)),  8, bouts_temp.onsets_loomaligned, '.', 'MarkerFaceAlpha', .2)
    deltas(idx_window, :) = avg_sm_window.(sprintf('avg_sm_%d', -current_window)) - bouts_temp.avg_sm;
    axis square
    apply_generic(gca, 'xlim', [0 20], 'ylim', [0 20])
    clim([0 50])

end
nexttile
scatter(bouts_temp.avg_sm, avg_sm_window.avg_sm_freeze, 12, '.')
axis square
apply_generic(gca, 'xlim', [0 20], 'ylim', [0 20])

colormap(cbrewer2('Spectral'))

%%
fh = figure('color','w','Position',[100, 100, 900, 400]);
tiledlayout(1, 5, 'TileSpacing', 'compact', 'Padding', 'loose')
for idx_mov_flies = 0:4
    nexttile
    mask = bouts_temp.moving_flies == idx_mov_flies;
    histogram(deltas(4, mask), -30:.5:30, 'Normalization', 'pdf', 'Orientation', 'horizontal')
    hold on
    apply_generic(gca, 'xlim', [0 1], 'ylim', [-30 30])
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