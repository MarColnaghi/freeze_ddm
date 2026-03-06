% Load the table with the computed log likelihoods 
model_list = {'dddim2'};
run_list   = {'run01_010326_1800'};
version = {'_v2'};

paths_analysis = path_generator('folder', fullfile('momentary_integration','alternative_test', model_list{1}, strcat(run_list{1}, version{1})));
col = cmapper('', 2);

ll = struct2table(importdata(fullfile(paths_analysis.results, 'll_table.mat')));
data = alternative_test_collect_data(model_list, run_list);

fh = figure('color', 'w','Position',[100, 100, 500, 500]);
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact')
nexttile

hold on
ylim([-2 2])
xlim([-50 height(ll) + 50])
ax = gca;

[sorted_deltall, idx_deltall_tv] = sort(ll.tv - ll.st);
idx_crossing = find(sorted_deltall > 0, 1);

fh = figure('color','w','Position',[100, 100, 500, 500]);
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact')
nexttile 
hold on

% 1. Define the limits for the fill area
min_val = -9.5; % Slightly lower than your axis limit
max_val = 2.5;  % Slightly higher than your axis limit

% 2. Fill the Upper Triangle (Model TV is better)
% Vertices: (min, min), (min, max), (max, max)
fill([min_val, min_val, max_val], [min_val, max_val, max_val], ...
    hex2rgb(col.timevarying_sm), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% 3. Fill the Lower Triangle (Model ST is better)
% Vertices: (min, min), (max, min), (max, max)
fill([min_val, max_val, max_val], [min_val, min_val, max_val], ...
    hex2rgb(col.stationary_sm), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% 4. Plot the reference line and scatter on top
plot([min_val max_val], [min_val max_val], 'k--', 'LineWidth', 1)
scatter(ll.st, ll.st_theory, 25, 1:length(ll.st), 'filled', 'MarkerFaceAlpha', 0.2)

% 5. Adjust labels and axis
axis square
xlim([min_val, max_val])
ylim([min_val, max_val])
xlabel('stationary LL')
ylabel('time-varying LL')
apply_generic(gca)


%%

% --- 1. Identify the outliers ---
% Calculate the absolute deviation from the diagonal
ll.deviation = abs(ll.st - ll.st_theory);

% Define a threshold (e.g., 0.5 units) or use the indices you found via inspection
% If you want to grab the specific points you saw (e.g. index 450 to 500):
problematic_idx = find(ll.deviation > 0.025); 

% Display the metadata for these specific bouts
disp('Bouts with significant deviation from theory:');
disp(data(1).table(problematic_idx, :));

% --- 2. Visualization with Annotated Deviations ---
fh = figure('color','w','Position',[100, 100, 600, 600]);
ax = axes('NextPlot', 'add', 'Box', 'off', 'TickDir', 'out');

% Background Triangles (Model comparison context)
min_val = -10; max_val = 3;
fill([min_val, min_val, max_val], [min_val, max_val, max_val], ...
    hex2rgb(col.timevarying_sm), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
fill([min_val, max_val, max_val], [min_val, min_val, max_val], ...
    hex2rgb(col.stationary_sm), 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Plot Diagonal
plot([min_val max_val], [min_val max_val], 'k:', 'LineWidth', 1.5)

% Scatter points - Color by Deviation to make them pop
scatter(ll.st, ll.st_theory, 30, ll.deviation, 'filled', 'MarkerFaceAlpha', 0.5);
colormap(ax, colorcet('R2')); % Linear highlight for the outliers
cb = colorbar; cb.Label.String = 'Abs Deviation from Theory';

% Highlight your "inspected" points specifically (optional)
plot(ll.st(problematic_idx), ll.st_theory(problematic_idx), 'ro', 'MarkerSize', 8);

% Labels
xlabel('Stationary LL (Computed)');
ylabel('Stationary LL (Theoretical)');
title('Verification of Stationary Model Likelihoods');
axis square; grid on;
xlim([min_val, max_val]); ylim([min_val, max_val]);
apply_generic(ax);