% Load the table with the computed log likelihoods
model_list = {'dddim2'};
run_list   = {'run01_010326_1800'};
version = {'_v3'};

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
min_val = 0; % Slightly lower than your axis limit
max_val = 2;  % Slightly higher than your axis limit

% 4. Plot the reference line and scatter on top
scatter(abs(ll.tv - ll.st), ll.js_div, 25, data.table.durations_s, 'filled', 'MarkerFaceAlpha', 0.4)
clim([0 5]);
colorcet('I1'); % Requires colorcet toolbox

% 5. Adjust labels and axis
xlim([min_val - 0.05, max_val + 0.05])
ylim([-.005 .205])
xlabel('$\log \!\big(\mathcal{L}_{\mathrm{tv}} / \mathcal{L}_{\mathrm{st}}\big)$', 'Interpreter','latex');
ylabel('js divergence')
apply_generic(gca)