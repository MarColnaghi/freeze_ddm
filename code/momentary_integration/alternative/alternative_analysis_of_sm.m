% --- 1. Setup and Data Loading ---
model_list = {'dddim2'};
run_list   = {'run01_010326_1800'};
version = {'_v8'};

paths_analysis = path_generator('folder', fullfile('momentary_integration','alternative_test', model_list{1}, strcat(run_list{1}, version{1})));

col = cmapper('', 2);

% Load Likelihoods and Data
ll = struct2table(importdata(fullfile(paths_analysis.results, 'll_table.mat')));
data = alternative_test_collect_data(model_list, run_list);

% Load Motion Cache (from your analyse_sm logic)
% Assuming data(1).sim_params contains the path to motion_cache
motion_cache = importdata(fullfile(paths_analysis.cache_path, 'motion_cache.mat'));

% --- 2. Configuration ---
align_to = 'offset'; % 'onset' or 'offset'
trace_len = numel(data.t);    % Frames to show
selection_size = 250; % How many top/bottom bouts to include in the average

% Remove the bouts that are censored?

freezes_ref = data(1).table;

% Remove censored freezes
is_censored = freezes_ref.durations_s > data(1).results.points.censoring;
freezes_ref = freezes_ref(~ is_censored,:);
ll = ll(~ is_censored,:);
data(1).table = data(1).table(~ is_censored,:);
data(1).extra_tv = data(1).extra_tv(~ is_censored,:);

n_bouts     = height(freezes_ref);
diff = ll.tv - ll.st;
[~, sort_idx] = sort(diff);

% Define the two groups
bottom_idx = sort_idx(1:selection_size);                 % Best for Stationary (Red)
top_idx    = sort_idx(end-selection_size+1:end);         % Best for Time-Varying (Blue)

% --- 3. Process Traces (Matrix Allocation) ---
bottom_mat = nan(selection_size, trace_len);
top_mat    = nan(selection_size, trace_len);

if ~isfield(data(1).results.estimates_mean, 'tndt_intercept')
    data(1).results.estimates_mean.tndt_intercept = 0;
end

if strcmpi(align_to, 'onset'), x = 0:trace_len-1; xt = 0:2/data.dt:3000; xtl = xt*data.dt; else, x = -trace_len+1:0; xt = -3000:2/data.dt:0; xtl = xt*data.dt; end

% Fill Matrices

for i = 1:selection_size
    % Extract Bottom Group

    bottom_mat(i, :) = get_cropped_signal_local(data, bottom_idx(i), trace_len, align_to);
    % Extract Top Group
    top_mat(i, :) = get_cropped_signal_local(data, top_idx(i), trace_len, align_to);
end

% --- 4. Main Plot ---
fh = figure('color','w','Position',[100, 100, 600, 450]);
hold on;

% Plot Individual Traces with high transparency (0.1)
plot(x, bottom_mat, 'Color', [hex2rgb(col.stationary_sm) 0.05], 'LineWidth', 0.5);
plot(x, top_mat, 'Color', [hex2rgb(col.timevarying_sm) 0.05], 'LineWidth', 0.5);

% Plot Grand Median Traces
plot(x, median(bottom_mat, 1, 'omitnan'), 'Color', col.stationary_sm, 'LineWidth', 2.5);
plot(x, median(top_mat, 1, 'omitnan'), 'Color', col.timevarying_sm, 'LineWidth', 2.5);

% Formatting
apply_generic(gca, 'xtick', xt); xticklabels(xtl);
y_min = -3; y_max = 3;
ylim([y_min, y_max]);
if strcmpi(align_to, 'onset'), xlim([0 trace_len]); else, xlim([-trace_len 0]); end

% Labels
if strcmpi(align_to,'onset')
    xlabel({'Time (seconds)', '[freeze start aligned]'});
    ylabel({'Soc. Motion Signal', '(start subtracted)'});
else
    xlabel({'Time (seconds)', '[freeze break aligned]'});
    ylabel({'Soc. Motion Signal', '(endpoint subtracted)'});
end

% Add a simple legend for the grand averages
% text(x(round(end*0.1)), 2.5, 'Best for Stationary', 'Color', col.stationary_sm, 'FontWeight', 'bold');
% text(x(round(end*0.1)), 2.1, 'Best for Time-Varying', 'Color', col.timevarying_sm, 'FontWeight', 'bold');

exporter(fh, paths_analysis, sprintf('grand_average_sm_%s.pdf', align_to));

% --- Helper for local execution ---
function cropped = get_cropped_signal_local(data, idx, trace_len, align_to)

signal = data.extra_tv(idx, :);
cropped = nan(1, trace_len);

dur_f = round(data.table.durations_s(idx) / data.dt);

if dur_f > numel(signal); return; end

seg = signal(1:dur_f);
n = min(numel(seg), trace_len);

if strcmpi(align_to,'onset')
    cropped(1:n) = seg(1:n);
    if ~isnan(cropped(1)); cropped = cropped - mean(cropped(1:2)); end
else
    cropped(end-n+1:end) = seg(end-n+1:end);
    if ~isnan(cropped(end)); cropped = cropped -  mean(cropped(end-1:end)); end
end
end