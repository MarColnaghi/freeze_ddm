% Fit bsl immobilities

clear all

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_freezes/bsl', 'bouts_id', id_code);
load(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'bsl', 'window', 'all');
points.censoring = 10.5;
points.truncation = min(bouts_proc.durations_s);

figure
hold on
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le');
histogram(bouts_proc.durations_s, 'BinEdges', points.truncation - 1/120 :1/20: points.censoring, 'Normalization', 'pdf')
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'bsl', 'window', 'all');
histogram(bouts_proc.durations_s, 'BinEdges', points.truncation - 1/120 :1/20: points.censoring, 'Normalization', 'pdf')
points.truncation = min(bouts_proc.durations_s);

model_results = run_fitting_newer(bouts_proc, points, 'exp0', paths, 'export', true, 'extra', []);
plot_fit('results', model_results, 'conditions', false)
