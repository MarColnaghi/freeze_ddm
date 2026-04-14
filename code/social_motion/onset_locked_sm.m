%% Freeze Onset aligned sm

bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le', 'nloom', 2:20, 'min_dur', 30);

window = [-300 630];
offsets = (window(1) : window(2));
total_looms = 20;
total_flies = size(motion_cache, 1);

n_moving_flies = nan(total_flies, 1);
sloom = nan(total_flies, 1);
sm_around_loom = nan(total_flies, total_looms, length(offsets));

for idx_fly = 1:total_flies
    fly_loom_x_sm = nan(total_looms, length(offsets));

    sm_fly = motion_cache(idx_fly);
    bouts_fly = bouts_proc(bouts_proc.fly == idx_fly, :);
    freeze_start = bouts_fly.onsets;
    freeze_end = bouts_fly.ends - bouts_fly.onsets - 1;

    idx_slice = freeze_start(:) + offsets;
    fly_loom_x_sm(bouts_fly.nloom, :) = sm_fly(idx_slice);

    try
        n_moving_flies(idx_fly) = unique(bouts_fly.moving_flies);
        sloom(idx_fly) = unique(bouts_fly.sloom);
    catch
    end
    
%     for i = 1:height(bouts_fly)
%         fly_loom_x_sm(bouts_fly.nloom(i), freeze_end(i):end) = nan;
%     end
    sm_around_loom(idx_fly, :, :) = fly_loom_x_sm;

end

col = cmapper( '', 30);
fh = figure('color','w','Position',[100, 100, 1400, 630]);
tiledlayout(2, 5, 'TileSpacing', 'compact')

graphvector.yaxis = [1 0 0 0 0 0 0 0 0 0];

sloom = sloom(~isnan(sloom));
n_moving_flies = n_moving_flies(~isnan(n_moving_flies));

for idx_sloom = unique(sloom)'

    for idx_moving_flies = unique(n_moving_flies)'

        nt = nexttile;
        hold on

        mask = n_moving_flies == idx_moving_flies & sloom == idx_sloom;
        sm_loom = squeeze(medi% 1. Create a VideoReader object
% Replace the filename with your specific .avi file
v = VideoReader('MH-1CS1NorpA3LC6Chrimson_5F-20LoomOpto100mW5Hz2min33cm-FH1_1-2021-04-09T11_38_31.avi');

% 2. Access the 'NumFrames' property
totalFrames = v.NumFrames;

% 3. Display the result
fprintf('The video has %d frames.\n', totalFrames);an(sm_around_loom(mask, :, :), 1, 'omitnan'));

        colororder(col.vars.ln(end - total_looms:end, :))
        plot(offsets, sm_loom')

        if idx_moving_flies == 0 && idx_sloom == 25
            xlabels = false;
            % xlabel('frames')
            text(-75, 0, '0', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Clipping', 'off', 'FontSize', 16)
            text(-75, 10, '10', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'Clipping', 'off', 'FontSize', 16)

        end


        plot([-70 -70], [0 10], 'k-', 'LineWidth', 2, 'Clipping', 'off');
        colororder(col.vars.ln(end - total_looms:end, :))

        apply_generic(gca, 'no_y', true, 'ylims', [-2 45], 'xlims', [-60 180], ...
            'no_xlabels', xlabels, 'xticks', [-60 -30 0 30 60 120 180], 'tick_length', 0.04, 'font_size', 18)

        plot([0 0], [-1 31], 'k-.', 'LineWidth', 1);
        text(1, 34, 'Freeze Onset', 'FontSize', 18)
        xtickangle(0);

        %         if idx_sloom == 25
        %             fill([0 0 30 30], [-3 -4 -4 -3], 'k', 'Clipping', 'off', 'EdgeColor', 'none')
        %             fill([thresholds.le_window_sl fliplr(thresholds.le_window_sl)], [-6 -6 -5 -5], 'b', 'Clipping', 'off', 'EdgeColor', 'none')
        %
        %         elseif idx_sloom == 50
        %             fill([0 0 30 30], [-3 -4 -4 -3], 'k', 'Clipping', 'off', 'EdgeColor', 'none')
        %             fill([thresholds.le_window_fl fliplr(thresholds.le_window_fl)], [-6 -6 -5 -5], 'b', 'Clipping', 'off', 'EdgeColor', 'none')
        %
        %         end

        if idx_sloom == 25
            xlabels = true;
            txtpos = nt.Position;
            txtpos(2) = txtpos(2) + 0.03;
            annotation('textbox', txtpos, 'String', sprintf('Moving Flies: %d', idx_moving_flies) ,'HorizontalAlignment', 'left', 'VerticalAlignment','top', ...
                'LineStyle', 'none', 'FontSize', 18);

        end

    end
end

exporter(fh, paths, 'freeze_onset_sm_std.pdf')



