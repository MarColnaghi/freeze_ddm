
clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'social_motion/clustering', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

% Here you can change the window of loom-evoked
%thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);

% Process the data: set thresholds for sm, fs and ln. Set minimum duration.
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);

%% Extract the pre-freezing social motion
sm_freeze = extract_sm_from_bouts(bouts_proc, 'type', 'onlyfreeze', 'output_type', 'mat');
inizio = 0;
fine = 240;
sm_freeze_full = extract_sm_from_bouts(bouts_proc, 'type', 'onsets', 'output_type', 'mat', 'window', [inizio fine]);

fh = figure('color','w','Position',[100, 100, 1800, 800]);
tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'tight')


for idx_loom_speed = 0:1

    mask = bouts_proc.ls == idx_loom_speed;
    bouts_sloom = bouts_proc(mask, :);
    sm_freeze_sloom = sm_freeze(mask, :);

    idx_loop = idx_loom_speed + 1;


    [~, sort_by_onset] = sort(bouts_sloom.onsets_loomaligned);

    nexttile(idx_loop)
    ax(idx_loop) = gca;
    hold on
    x_axis = 1:size(sm_freeze_sloom, 1);
    y_axis = inizio:size(sm_freeze_sloom,2);

    h = imagesc(x_axis, y_axis, sm_freeze_sloom(sort_by_onset, :)', [0 0.4]);

    % set(h, 'AlphaData', ~isnan(sm_freeze(sort_order, :)));

    scatter(1:height(bouts_sloom), bouts_sloom.durations(sort_by_onset),  2, '_', 'k')
    apply_generic(gca, 'ylim', [0 630], 'no_x', true, 'xlim', [-100 size(sm_freeze_sloom, 1) + 100], 'ytick', 0:120:600)
    colormap(cbrewer2('Reds',[]));
end

%cb = colorbar(ax, 'Location', 'southoutside', 'FontSize', 18, 'LineWidth', 2);
%cb.Label.String = 'Social Motion';

% for idx_cluster = 1:n_clusters
%     fill([ax(2).XLim(1)-2,ax(2).XLim(1)-2,ax(2).XLim(1)-17,ax(2).XLim(1)-17], [boundaries(idx_cluster) boundaries(idx_cluster + 1) boundaries(idx_cluster + 1)   boundaries(idx_cluster)], ...
%         col.pca(idx_cluster,:), 'EdgeColor','none', 'Clipping', 'off');
%     text(mean([ax(2).XLim(1)-2,ax(2).XLim(1)-17]), mean([boundaries(idx_cluster); boundaries(idx_cluster + 1)]), num2str(idx_cluster),...
%         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
%     hold on
% 
% end

% xlabel('Time');
% 
% % Add horizontal lines to separate clusters
% hold on;
% for idx_cluster = 1:length(boundaries)
%     yline(boundaries(idx_cluster), 'k--', 'LineWidth', 1.1);
% end
% 
% nexttile(1)
% ax(1) = gca;
% hold on
% for idx_cluster = 1:n_clusters
%     plot(inizio:fine, repr(idx_cluster, :), 'LineWidth', 1.5, 'Color', col.pca(idx_cluster, :));
% end
% apply_generic(gca, 'ylim', [0 1.2],'ytick', [0 1.2], 'xlim', [0 630], 'no_x', true, 'font_size', 18);
% ylabel({'Mean', 'Social Motion'})

linkaxes([ax(:)], 'xy')