% Fit bsl immobilities

clear all

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 3; threshold_mob = 3; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_freezes/bsl', 'bouts_id', id_code);
load(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'loom', 'window', 'all');

figure
hold on
histogram(bouts_proc.durations_s, 'BinEdges', 1/120:1/60:10)
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'bsl', 'window', 'all');
histogram(bouts_proc.durations_s, 'BinEdges', 1/120:1/60:10)
