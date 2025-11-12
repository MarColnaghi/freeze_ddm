include_old = true;

paths = path_generator('folder', 'scoring');
pixel_cache = importdata(fullfile(paths.cache_path, 'pixel_cache.mat'));
loom_cache = importdata(fullfile(paths.cache_path, 'loom_cache.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));
thresholds = define_thresholds;
fly_id = 862;
fly_pc = ~pixel_cache(fly_id);
fly_loom = loom_cache(fly_id);
loom_frames = find(diff(fly_loom) == 1) + 1;

fly_sm = motion_cache(fly_id);
fh = figure('color','w','Position',[100, 100, 1600, 500]);
tiledlayout(4, 1)
nexttile
imagesc(fly_pc', [0 6])
ax(1) = gca;
apply_effects(ax(1));

if include_old
    X = importdata(fullfile('/Users/marcocolnaghi/PhD/social_fly_ddm/datasets', 'X1.mat'));
    xs = filter_x(X);
    xs_fly = xs(xs.fly == fly_id, :);
    old_fly_pc = zeros(1, length(fly_sm));
    for idx_freezes = 1:height(xs_fly)
        starts = xs_fly.freeze_start(idx_freezes); ends = xs_fly.freeze_start(idx_freezes) + round(xs_fly.freeze_dur(idx_freezes) * 60); 
        old_fly_pc(starts:ends) = 2;
    end
    nexttile
    hold on
    ax(4) = gca;
    apply_effects(ax(4));
    imagesc(old_fly_pc, [0 2])

end

id_code = 'imm2_mob2_pc4';
types = {'fill_in_combinations_mobfirst', 'fill_in_combinations_imfirst'};


for idx_bout_solution = 1:2

    bouts = importdata(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/datasets', types{idx_bout_solution}, id_code,  'bouts.mat'));

    thresholds = define_thresholds;
    %thresholds.le_window_fl = [5 40];
    %thresholds.le_window_sl = [15 50];

    bouts = bouts_formatting(bouts, thresholds);
    bouts_fly = bouts(bouts.fly == fly_id, :);
    bouts_fly = sortrows(bouts_fly, {'onsets'});
    dur = bouts_fly.durations;
    dur(dur <= 0 | isnan(dur)) = 0;
    base_val = double(bouts_fly.type);
    base_val(bouts_fly.le == true & bouts_fly.type == true & bouts_fly.frozen_start == false & bouts_fly.durations >= 30) = 2;

    ts_all = repelem(base_val, dur);

    nexttile
    hold on
    ax(1 + idx_bout_solution) = gca;
    apply_effects(ax(1 + idx_bout_solution));
    imagesc(ts_all', [0 2])

end

plot([loom_frames loom_frames], [7 -2], 'k', 'Clipping', 'off', 'LineStyle', '--');
if bouts_fly.sloom(1) == 25
    fts = loom_frames + thresholds.le_window_sl;
elseif bouts_fly.sloom(1) == 50
    fts = loom_frames + thresholds.le_window_fl;
end
plot([fts(:, 1) fts(:, 1)] , [7 -2], 'w', 'Clipping', 'off', 'LineStyle', '--');
plot([fts(:, 2) fts(:, 2)] , [7 -2], 'w', 'Clipping', 'off', 'LineStyle', '--');

linkaxes(ax)
ylim([0.5 1.5])
xlim([0 length(fly_pc)])

function apply_effects(ax)
axis off
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
end

function Xs = filter_x(X)
Nloom = mod(0:height(X)-1,20)'+1;
X.Nloom = Nloom;
X_short = X(X.hit_looms==1,:);
X_short.freeze_dur = X_short.freeze_dur/60;
X_short.loom_interval = X_short.loom_interval/60;
X_short.freeze_latency = X_short.freeze_latency/60;
X_short.avg_sm_freeze = X_short.avg_sm_freeze/10;
X_short.avg_fs_pre = X_short.avg_fs_pre/10;
X_short.avg_fs_1s_prefreeze = X_short.avg_fs_1s_prefreeze/10;
X_short.loom_speed = X_short.loom_speed/25;
X_short.Nloom = X_short.Nloom/10;
X_short = X_short(X_short.freeze_dur >= 0.5,:);
X_short = X_short(X_short.freeze_latency<=1,:);
Xs = X_short(X_short.avg_sm_freeze < 1.25 &...
    X_short.avg_fs_1s_prefreeze < 2,:);
end