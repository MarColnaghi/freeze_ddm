

clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'collision_analysis', 'bouts_id', id_code, 'imfirst', false);
bouts = importdata(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));
sumangle_cache = importdata(fullfile(paths.cache_path, 'sumangle_cache.mat'));
sumspeed_cache = importdata(fullfile(paths.cache_path, 'sumspeed_cache.mat'));
mindist_cache = importdata(fullfile(paths.cache_path, 'mindist_cache.mat'));

surangle_cache = importdata(fullfile(paths.cache_path, 'surangle_cache.mat'));
surspeed_cache = importdata(fullfile(paths.cache_path, 'surspeed_cache.mat'));
surdistance_cache = importdata(fullfile(paths.cache_path, 'surdistance_cache.mat'));

% Here you can change the window of loom-evoked
thresholds = define_thresholds('le_window', struct('le_window_sl', [5 55], 'le_window_fl', [5 55]));
% thresholds = define_thresholds;

bouts = bouts_formatting(bouts, thresholds);

%%

idx_fly = 13;

% --- 1. Data Extraction ---
sm     = motion_cache(idx_fly);       % Note: currently unused in plots
ss_1   = surspeed_cache(idx_fly);     % Note: currently unused in plots
sa_1   = surangle_cache(idx_fly);
dist_1 = surdistance_cache(idx_fly);

idx_surr_fly = 1;

% --- 2. Timeseries Plot (Distance & Angle) ---
fig1 = figure('Color', 'w', 'Position', [100, 100, 1200, 400], 'Name', 'Timeseries');
tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'tight');
ax1 = nexttile;
hold(ax1, 'on');

% Right Axis: Angle
yyaxis(ax1, 'right');
plot(ax1, sa_1(:, idx_surr_fly), 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
ylabel(ax1, 'Angle');

% Left Axis: Distance
yyaxis(ax1, 'left');
plot(ax1, dist_1(:, idx_surr_fly), 'LineWidth', 1.5, 'Color', [0 0.4470 0.7410]);
ylabel(ax1, 'Distance');
xlabel(ax1, 'Time (Frames)');

% --- 3. Scatter Plot (Angle vs. Distance) ---
fig2 = figure('Color', 'w', 'Position', [100, 100, 900, 500]); % Slightly enlarged to 500x500 for the inset
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'tight');
ax2 = nexttile;
hold(ax2, 'on');

% Find out how many surrounding flies we have (number of columns)
num_surr_flies = size(sa_1, 2); 

% Preallocate a cell array for the legend entries
legend_labels = cell(1, num_surr_flies);

% Loop through each surrounding fly column
for j = 1:num_surr_flies
    % Scatter plot for surrounding fly 'j'
    scatter(ax2, sa_1(:, j), dist_1(:, j), 8, 'filled', 'MarkerFaceAlpha', 0.25);
    
    % Create a label for the legend (e.g., "Fly 1")
    legend_labels{j} = sprintf('Fly %d', j);
end

ylabel(ax2, 'Distance');
xlabel(ax2, 'Angular Size');

% Placed legend explicitly in the northwest so the inset can take the top right
legend(ax2, legend_labels, 'Location', 'northwest', 'FontSize', 18, 'box', 'off');

% Apply your custom formatting to the main scatter plot
apply_generic(gca, 'font_size', 20, 'tick_length', 0.015, 'line_width', 2, 'xlim', [0 pi], 'ylim', [0 900]);

% xlims = [0.415, 0.417];
% ylims = [118, 119];
% 
% % Define the X and Y coordinates for the 4 corners of the box
% x_box = [xlims(1), xlims(2), xlims(2), xlims(1)];
% y_box = [ylims(1), ylims(1), ylims(2), ylims(2)];
% 
% % Draw the box (FaceAlpha = 0 makes the fill transparent)
% fill(ax2, x_box, y_box, 'k', 'FaceAlpha', 0, 'EdgeColor', 'k', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% % --- 4. Inset Axis (Zoom) ---
% % Define the inset position [left bottom width height] in normalized coordinates
% ax_inset = axes('Position', [0.65 0.65 0.25 0.25]); 
% hold(ax_inset, 'on');
% box(ax_inset, 'on'); % Add a border around the inset
% 
% % Re-plot the exact same data into the inset
% for j = 1:num_surr_flies
%     scatter(ax_inset, sa_1(:, j), dist_1(:, j), 8, 'filled', 'MarkerFaceAlpha', 0.25);
% end
% 
% xlim(ax_inset, xlims);  % Example: zooming in on angle between -0.5 and 0.5
% ylim(ax_inset, ylims);       % Example: zooming in on distance between 0 and 5
% 
% % Format the inset (using your custom function but with slightly smaller font/lines)
% apply_generic(ax_inset, 'font_size', 14, 'tick_length', 0.015, 'line_width', 2, 'box', 'on');
% 
% % Determine the total number of focal flies in your cache
% % (Assuming surangle_cache is a cell array or array of matrices)
% num_focal_flies = length(surangle_cache); 
% 
% % Zoom region for all plots
% xlims = [0.415, 0.417];
% ylims = [118, 119];
% x_box = [xlims(1), xlims(2), xlims(2), xlims(1)];
% y_box = [ylims(1), ylims(1), ylims(2), ylims(2)];
% 
% hold(ax_inset, 'off');
% box(ax_inset, 'off'); % Add a border around the inset

figure(fig2)
nexttile
hold on

% Preallocate arrays to store all points
num_flies = size(surangle_cache, 1);
samples_per_fly = 50;

all_sa = NaN(num_flies * samples_per_fly, 1);
all_dist = NaN(num_flies * samples_per_fly, 1);
all_id = NaN(num_flies * samples_per_fly, 1);

for k = 1:num_flies
    sa_k   = surangle_cache(k);
    dist_k = surdistance_cache(k);
    
    % Sample 100 random row indices
    idx_samples = randi(size(sa_k, 1), samples_per_fly, 1);
    
    % Calculate the indices for our preallocated array
    idx_start = (k - 1) * samples_per_fly + 1;
    idx_end   = k * samples_per_fly;
    
    % Store the data
    all_sa(idx_start:idx_end)   = sa_k(idx_samples, 1);
    all_dist(idx_start:idx_end) = dist_k(idx_samples, 1);
    all_id(idx_start:idx_end)   = k; % Tag these rows with the fly ID
end

% Plot EVERYTHING in a single command
scatter(all_sa, all_dist, 5, all_id, 'filled', 'MarkerFaceAlpha', 1);

ylabel('Distance');
xlabel('Angular Size');
apply_generic(gca, 'font_size', 20, 'tick_length', 0.015, 'line_width', 2);
exporter(fig2, paths, 'distance_vs_angle.png')
% exporter(fig2, paths, 'distance_vs_angle.pdf')

%% Find functional relationship that links distance to angle

% Define the fitting equation based on the paper
% 'a' is the unknown size, 'd' is distance, 'theta' is angle
ft = fittype('2 * atan(a / (2 * d))', 'independent', 'd', 'dependent', 'theta');

% Preallocate arrays to store the results for each fly
num_flies = size(surangle_cache, 1);
estimated_sizes = NaN(num_flies, 1);
r_squared       = NaN(num_flies, 1);

% Generate a color palette for the lines to match the scatter plot
colors = lines(num_flies);

% Set up the figure
figure('Color', 'w', 'Position', [100, 100, 600, 300])
hold on;

% 1. Optional: Plot the background scatter of ALL data in light gray for context
% (Assuming 'all_sa' and 'all_dist' are still in your workspace from earlier)
scatter(all_sa, all_dist, 5, [0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha', 0.1);

% 2. Loop through each fly to fit and plot
for k = 1:num_flies
    % Extract data for focal fly 'k', looking at surrounding fly 1
    sa_k       = surangle_cache(k);
    dist_k     = surdistance_cache(k);
    sa_surr1   = sa_k(:, 1);
    dist_surr1 = dist_k(:, 1);
    
    % Clean NaNs for this specific fly
    valid_idx  = ~isnan(sa_surr1) & ~isnan(dist_surr1);
    clean_sa   = sa_surr1(valid_idx);
    clean_dist = dist_surr1(valid_idx);
    
    % Only attempt to fit if there are enough valid data points (e.g., > 10)
    if length(clean_sa) > 10
        % Run the fit with a starting guess for size 'a' (e.g., 2)
        [fit_result, gof] = fit(clean_dist, clean_sa, ft, 'StartPoint', 2);
        
        % Store the results
        estimated_sizes(k) = fit_result.a;
        r_squared(k)       = gof.rsquare;
        
        % Plot the fitted line for this fly
        dist_curve   = linspace(min(clean_dist), max(clean_dist), 100);
        theta_fitted = fit_result(dist_curve);
        
        plot(theta_fitted, dist_curve, '-', 'Color', colors(k, :), 'LineWidth', 1.5);
    else
        warning('Fly %d does not have enough valid data points to fit.', k);
    end
end

% Format the plot
ylabel('Distance');
xlabel('Angle');
title('Individual Fits per Fly');
apply_generic(gca, 'font_size', 20, 'tick_length', 0.015, 'line_width', 2);

% --- Display Summary Statistics ---
fprintf('\n--- Fit Summary Across %d Flies ---\n', num_flies);
fprintf('Mean Estimated Size: %.3f\n', mean(estimated_sizes, 'omitnan'));
fprintf('Median Estimated Size: %.3f\n', median(estimated_sizes, 'omitnan'));
fprintf('Standard Deviation of Size: %.3f\n', std(estimated_sizes, 'omitnan'));
fprintf('Mean R-squared: %.3f\n', mean(r_squared, 'omitnan'));


% Set up the figure
figure('Color', 'w', 'Position', [100, 100, 600, 600])
hold on; 
histogram(estimated_sizes, 45:0.05:55)
apply_generic(gca, 'ylim', [-10 210], 'xlim', [45 55])