function fh = plot_sta_sm(aligned_mat, varargin)

opt = inputParser;
addParameter(opt, 'direction', 'offset');
addParameter(opt, 'center', 'median');
addParameter(opt, 'exporting', false);
addParameter(opt, 'clim', [0 15]);

parse(opt, varargin{:});

fh = figure('color','w','Position',[100, 100, 600, 650]);
tl = tiledlayout(4, 1, 'TileSpacing', 'tight', 'Padding', 'compact');

camera_fps = 60;

if strcmp(opt.Results.direction, 'offset')
    xvec = -10:2:0;
    lines = size(aligned_mat, 2) + xvec * camera_fps;
elseif strcmp(opt.Results.direction, 'onset')
    xvec = 10:-2:0;
    lines = 1 + xvec * camera_fps;
end

x = 1:size(aligned_mat, 2);

nexttile
hold on

% Calculate important statistics
if strcmp(opt.Results.center, 'median')
    mean_array = median(aligned_mat, 1, 'omitnan');
elseif strcmp(opt.Results.center, 'mean')
    mean_array = mean(aligned_mat, 1, 'omitnan');
end

sem_array = std(aligned_mat, 'omitnan') ./ sqrt(sum(~isnan(aligned_mat)));

fill([x fliplr(x)], ...
    [mean_array + sem_array fliplr(mean_array - sem_array)], ...
    'k', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
plot(mean_array, 'k', 'LineWidth', 3)
ax(1) = gca;
ax(1).XAxis.Visible = "off";
ylabel({'Social','Motion'})
apply_generic(ax(1), 20)
 
xline(lines, '--', 'LineWidth', 0.5)
ylim(opt.Results.clim)

nexttile(2,[3 1])
imgsc = imagesc(aligned_mat, opt.Results.clim);
set(imgsc, 'AlphaData', ~isnan(aligned_mat));
yticks([])

%ylim([0 max([length(early_sm), length(late_sm)])])

ax(2) = gca;
apply_generic(ax(2), 20);

set(gca,'YAxisLocation', 'right')
set(gca, 'XAxisLocation', 'top')

xlabel(sprintf('Time from Freeze %s (s)', opt.Results.direction))
ylabel('Freezes')
colorcet('L08')
linkaxes([ax(2), ax(1)], 'x')

if strcmp(opt.Results.direction, 'offset')
    xl = [size(aligned_mat, 2) + xvec(1) * camera_fps, size(aligned_mat, 2)];
    xlim(ax(1), xl);
    xlim(ax(2), xl);    
    set(ax(1),'YAxisLocation', 'right')
    set(ax(2),'YAxisLocation', 'right')
    xticks(sort(lines))
    xticklabels(sort(xvec))

elseif strcmp(opt.Results.direction, 'onset')
    xl = [0, abs(xvec(1)) * camera_fps];
    xlim(ax(1), xl);
    xlim(ax(2), xl);
    set(ax(1),'YAxisLocation', 'left')
    set(ax(2),'YAxisLocation', 'left')
    xticks(sort(lines))
    xticklabels(sort(xvec))

end

