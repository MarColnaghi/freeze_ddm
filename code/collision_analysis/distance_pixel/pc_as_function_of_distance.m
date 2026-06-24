clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'collision_analysis', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
mindist_cache = importdata(fullfile(paths.cache_path, 'mindist_cache.mat'));
pc_cache = importdata(fullfile(paths.cache_path, 'pixel_cache.mat'));

%%
% Define periods to analyze
period_names = {'Baseline', 'Loom'};
period_slices = {1:18000, 18001:36000}; 

% Initialize a structure to hold results
analysis_output = struct();

for p = 1:length(period_names)
    name = period_names{p};
    slice = period_slices{p};
    
    % Data Extraction
    dist_ts = [];
    pc_ts = [];
    for idx_flies = 1:size(pc_cache, 1)

        temp_md = mindist_cache(idx_flies);
        temp_pc = pc_cache(idx_flies);

        dist_ts = [dist_ts; temp_md(slice)'];
        pc_ts = [pc_ts; temp_pc(slice)'];
    end
    
    % Run the math
    fprintf('Computing threshold for %s...\n', name);
    params = struct('edges', 0:2.5:1000);
    res.dist_raw = dist_ts; % Save all distance points for this period
    res = compute_distance_threshold(pc_ts, dist_ts, params);

    % SAVE THE DISTRIBUTION DATA
    res.dist_raw = dist_ts; % Save all distance points for this period
    analysis_output.(name) = res;
    fprintf('Done computing threshold for %s...\n', name);

end

%%

% Define styles for each period
styles = struct(...
    'Baseline', struct('color', [0.2, 0.2, 0.8], 'label', 'Pre-Stimulus'), ...
    'Loom',     struct('color', [0.8, 0.2, 0.2], 'label', 'Looming Period')...
);

fn = fieldnames(analysis_output);
h_lines = []; % For the legend

fh = figure('color','w','Position',[100, 100, 500, 800]);
tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'tight')

nexttile
hold on;
ax(1) = gca;

for i = 1:length(fn)
    res = analysis_output.(fn{i});
    cfg = styles.(fn{i});
    histogram(res.dist_raw, params.edges, 'Normalization', 'pdf', 'FaceColor', cfg.color, 'EdgeColor', 'none')
end

apply_generic(gca, 'font_size', 24, 'tick_length', 0.015, 'line_width', 2, 'xlim', [0 900]);
xline(60, '--', 'Color', 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');

nexttile(2, [2,1])
hold on;
ax(2) = gca;

for i = 1:length(fn)
    res = analysis_output.(fn{i});
    cfg = styles.(fn{i});
    
    % Plot raw binned data (semi-transparent)
    line(res.binCenters, res.P_bin, 'LineStyle','none', 'Marker','.', 'MarkerSize', 8, 'MarkerEdgeColor', cfg.color, ...
         'HandleVisibility', 'off');
    
    % Plot smooth curve
    h_lines(i) = plot(res.D_fit, res.P_smooth, 'Color', cfg.color, 'LineWidth', 1.5);
    
    % Plot threshold line
    xline(60, '--', 'Color', cfg.color, 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

% Polish the look

xlabel('Distance');
ylabel('Pixel Change');
legend(h_lines, {styles.Baseline.label, styles.Loom.label}, 'Location', 'northeast', 'box', 'off', 'FontSize', 20);
apply_generic(gca, 'font_size', 24, 'tick_length', 0.015, 'line_width', 2, 'xlim', [0 900]);
linkaxes(ax(:), 'x')


exporter(fh, paths, 'binned_distance.pdf')
% exporter(fh, paths, 'binned_distance.png')

%% Here we can do a Single-Fly plot 

params = struct('edges', 0:1:1000);
fh = figure('color','w','Position',[100,100,500,800]);
apply_generic(gca,'font_size',24,'tick_length',0.015,'line_width',2,'xlim',[0 900],'ylim',[0 900]);

% --- 1. Initialize Video ---
v = VideoWriter(fullfile(paths.fig, 'fly_distance_highlight.mp4'),'MPEG-4');
v.FrameRate = 10;
open(v);
hold on

% --- Pre-create the two reusable objects (created ONCE) ---
sh_hist = line(nan, nan, 'LineStyle','none', 'Marker','.', 'MarkerSize',8, ...
    'Color',[0 0 0], 'HandleVisibility','off');
sh_curr = line(nan, nan, 'LineStyle','none', 'Marker','.', 'MarkerSize',12, ...
    'Color',[1 0 0], 'HandleVisibility','off');

xlabel('Distance');
ylabel('Pixel Change');

% Buffers for the accumulated (black) points
nFlies   = size(pc_cache, 1);
hist_x   = [];
hist_y   = [];
prev_x   = [];   % last fly's points, to be folded into history
prev_y   = [];

for idx_flies = 1:nFlies
    temp_md = mindist_cache(idx_flies);
    temp_pc = pc_cache(idx_flies);
    dist_ts = temp_md(:);
    pc_ts   = temp_pc(:);

    % Run the math
    res_single = compute_distance_threshold(pc_ts, dist_ts, params);
    analysis_output_single.('fly_id')(idx_flies) = res_single;

    % --- Fold the previous current-fly into the black history cloud ---
    if idx_flies >= 2
        hist_x = [hist_x; prev_x];
        hist_y = [hist_y; prev_y];
        set(sh_hist, 'XData', hist_x, 'YData', hist_y);
    end

    % --- Update the current (red) fly in place ---
    prev_x = res_single.binCenters(:);
    prev_y = res_single.P_bin(:);
    set(sh_curr, 'XData', prev_x, 'YData', prev_y);

    % --- Capture the frame ---
    drawnow;
    frame = getframe(fh);
    writeVideo(v, frame);
end


% --- Clean up ---
close(v);
hold off;

%%





figure
only_blind_flies = unique(bouts.fly(bouts.moving_flies == 4, :));
only_dist = dist_ts(only_blind_flies, :);
only_pc = pc_ts(only_blind_flies, :);

scatter(only_dist(:), only_pc(:), 2, '.', 'k')
