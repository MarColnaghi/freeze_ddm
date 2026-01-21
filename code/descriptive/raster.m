clearvars

thresholds.fill_in_imm = 2;
thresholds.fill_in_mob = 2;
thresholds.pc = 4;
thresholds.imfirst = false;

id_code = sprintf('imm%d_mob%d_pc%d', thresholds.fill_in_imm, thresholds.fill_in_mob, thresholds.pc);
paths = path_generator('folder', 'descriptive/raster', 'bouts_id', id_code, 'imfirst', false);
col = cmapper();
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
sm_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));
pc_cache = importdata(fullfile(paths.cache_path, 'pixel_cache.mat'));
fs_cache = importdata(fullfile(paths.cache_path, 'speed_cache.mat'));

% Select Condition

n_looms = 20;
loom_speed = 50;
n_moving_flies = 3;

% Minimum duration freeze
min_dur_freeze = 0.5;

% Select only relevant flies

bouts_loom =  data_parser_new(bouts, 'period', 'loom');
loom_ts = nan(n_looms, 1);

% Identify Loom Timestamps
for idx_loom = 1:n_looms
   
    loom_ts(idx_loom) = round(mean(bouts_loom.loom_ts(bouts_loom.nloom == idx_loom)));

end

flies = unique(bouts.fly(bouts.moving_flies == n_moving_flies & bouts.sloom == loom_speed));
n_flies = size(flies, 1);

% Preallocate the arrays
sm_flies = nan(n_flies, 36000);
pc_flies = nan(n_flies, 36000);
fs_flies = nan(n_flies, 36000);

% Prepare the arrays
for idx_fly = 1:size(flies, 1)
    
    fly_number = unique(bouts.fly(bouts.fly == idx_fly));
    
    sm_flies(idx_fly, :) = sm_cache(fly_number);   
    pc_flies(idx_fly, :) = pc_cache(fly_number);
    fs_flies(idx_fly, :) = fs_cache(fly_number);

end

%

thresholds.fill_in_imm = 30;
thresholds.fill_in_mob = 2;
thresholds.pc = 4;

imm_frames = pc_flies < thresholds.pc;
imm_frames = cell2mat(cellfun(@(row) ...
    processRow(row, thresholds), ...
    num2cell(imm_frames, 2), ...
    'UniformOutput', false));


fh = figure('color', 'w', 'Position', [100 200 800 700]);
tl = tiledlayout(1, 1, 'TileSpacing', 'tight', 'Padding', 'compact');


nexttile
ax2 = gca;
total_freezetime = sum(imm_frames, 2);

[sorted_by_total_freezetime, idx_by_total_freezetime] = sort(total_freezetime);

fs_flies(imm_frames == 1) = nan;
fs_flies = fs_flies(idx_by_total_freezetime, :);

sp_imgsc = imagesc(ax2, fs_flies, [0, 10]);
colorcet('L01')
set(sp_imgsc, 'AlphaData', ~isnan(fs_flies));


%
% 
% nexttile
% hold on
% ax1 = gca;
% ax1.YColor = 'k';
% ax1.XColor = 'k';
% ax1.FontSize = 32;
% yl = ylabel(ax1, {'Fraction', 'Freezing'}, 'FontSize', 32);
% xl = xlabel(ax1, 'Time (min)', 'FontSize', 32);
% set(ax1, 'Xtick', [0, loom_ts(1), length(imm_frames)], 'XTicklabels', [0, 5, 10], 'YTick', [0, 1]);
% set(ax1,'box','off');
% ax1.XAxis.Visible = 'off';
% set(ax1,'linewidth', 3, 'TickDir','out'); % The only other option is 'in'p
% 

function row = processRow(row, thresholds)
    if thresholds.imfirst
        row = bwareaopen(row, thresholds.fill_in_imm);
        row = ~bwareaopen(~row, thresholds.fill_in_mob);
    else
        row = ~bwareaopen(~row, thresholds.fill_in_mob);
        row = bwareaopen(row, thresholds.fill_in_imm);
    end
end
