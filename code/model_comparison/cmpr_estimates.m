clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
model_str = 'dddm2';
model_func = str2func(strcat('model_',model_str));
model_struct = model_func();
col = cmapper();
col.Set3 = cbrewer2('Set3', 5);

markers = {'o','+','*','s','d','v','>','h'};

mkrsize = 220;

model_list = {'dddm2', 'ded2'};
run_list   = {'run21', 'run15'};
label_list = {'Integration', 'Extrema Detection'};
color_list = {col.timevarying_sm, col.extremadetection};

fh = figure('color','w', 'Position', [100 100 1200 400]);

[fh] = base_for_estimates('model', model_struct, 'ylimits', [-1 5], 'fh', fh);
paths = path_generator('folder', 'fitting_freezes/le');

for idx_model = 1:length(model_list)

    results_folder = fullfile(paths.results, model_list{idx_model}, run_list{idx_model});

    model_results = importdata(fullfile(results_folder, sprintf('fit_results_%s.mat', model_list{idx_model})));

    plot_estimates('base', fh, 'results', model_results, 'marker', markers{idx_model}, 'label', label_list{idx_model}, 'ylimits', [-1 5], 'colors', color_list{idx_model})
end

legend('Location', 'eastoutside', 'box', 'off', 'FontSize', 20)
paths = path_generator('folder', 'model_comparison');
exporter(fh, paths, 'comparison.pdf')