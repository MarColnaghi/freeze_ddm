% Double DDM test

clear all

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_freezes/le', 'bouts_id', id_code, 'imfirst', false);
load(fullfile(paths.dataset, 'bouts.mat'));
motion_cache = importdata(fullfile(paths.cache_path, 'motion_cache.mat'));

bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'le');
points.censoring = 10.5;
points.truncation = 0.3;

%%  Now we added our vector column to the bouts table

model_results = run_fitting_newer(bouts_proc, points, 'dddm2', paths, 'export', true, 'extra', []);
model_results.estimates_mean  
plot_fit('results', model_results, 'conditions', false)
plot_fit('results', model_results, 'conditions', true)
plot_estimates('results', model_results, 'export', true, 'paths', paths)
