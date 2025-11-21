clearvars
paths = path_generator('folder', '/social_motion');
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));
pixel_cache = importdata(fullfile(paths.cache_path, 'pixel_cache.mat'));
loom_cache = importdata(fullfile(paths.cache_path, 'loom_cache.mat'));
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));

window = [-300 300];
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

%% Plot Time Locked Social Motion

col = cmapper( '', 30);
fh = figure('color','w','Position',[100, 100, 1400, 400]);
tiledlayout(2, 5)

graphvector.yaxis = [1 0 0 0 0 1 0 0 0 0];

for idx_sloom = unique(sloom)'

    for idx_moving_flies = unique(n_moving_flies)'
        nexttile
        mask = n_moving_flies == idx_moving_flies & sloom == idx_sloom;
        sm_loom = squeeze(mean(sm_around_loom(mask, :, :), 1));

        plot(offsets, sm_loom')
        colororder(col.vars.ln(end - total_looms:end, :))
        apply_generic(gca, 'no_y', true, 'ylims', [-2 45])
        xline(0, 'LineWidth', 2); xtickangle(0);
        xticks(sort([window, 0]))
    end
end

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
