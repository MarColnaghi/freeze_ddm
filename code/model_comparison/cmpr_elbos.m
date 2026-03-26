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

model_list = {'dddm2', 'dddim2', 'ded2','ded2'};
run_list   = {'run14_260325', 'run02_260325', 'run05_260325', 'run06_260325'};
label_list = {'ddm (with ndt)', 'ddm (no ndt)', 'extrema detection', 'extrema detection (constant evidence)'};
color_list = {col.timevarying_sm, col.extremadetection};
color_list = {'#20DA99', '#5921E6', col.extremadetection, col.extremadetection};

fh = figure('color','w', 'Position', [100 100 700 500]);

paths = path_generator('folder', 'fitting_freezes/le_new');
elbo = nan(length(model_list), 2);

for idx_model = 1:length(model_list)

    results_folder = fullfile(paths.results, model_list{idx_model}, run_list{idx_model});

    model_results = importdata(fullfile(results_folder, 'model_results.mat'));

    elbo(idx_model, : ) = model_results.elbo;

end

b = bar(1:length(model_list), elbo(:, 1), 'FaceColor', 'flat', 'EdgeColor', 'none');

% Assign colors from color_list to each individual bar
for k = 1:length(b.CData)
    b.CData(k,:) = validatecolor(color_list{k});
end

paths = path_generator('folder', 'model_comparison');
ax = gca;
apply_generic(ax, 'ylim', [-12000 -10000])
xticklabels(label_list);
exporter(fh, paths, 'elbo_wwo_ndt.pdf')

