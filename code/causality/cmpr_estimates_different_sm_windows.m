clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
model_str = 'dddm2';
model_func = str2func(strcat('model_',model_str));
model_struct = model_func();
col.Set3 = cbrewer2('Set3', 5);

windows.anchor = 'freeze_onset';
windows.reference = 'fixed_length';
windows.length = '30';
windows.select_windows = {'-300' ,'-60', '0' ,'10' ,'30' ,'60'};
markers = {'o','+','*','s','d','v','>','h'};

fh = figure('Position', [100 100 920 560], 'Color', 'w');
t = tiledlayout(3, 1, 'TileSpacing', 'loose', 'Padding', 'loose');
mkrsize = 220;

model_acausal = importdata('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/fitting_freezes/le/dddm2/run11/fit_results_dddm2.mat');
elbo_acausal = model_acausal.elbo;

parent_folder = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/causality/fitting_windows', model_str, windows.reference, 'run01_size30_anchor-freeze_onset');
d = dir(fullfile(parent_folder));
d = d(~ismember({d.name}, {'.','..', '.DS_Store'}));

get_window_size = nan(size(d, 1), 1);

[fh] = base_for_estimates('model', model_struct);
[~, idx] = ismember(windows.select_windows, {d.name});

i = 0;
for idx_windows = idx
    i = i + 1;
    get_window_size(idx_windows) = str2double(d(idx_windows).name);
    model_results = importdata(fullfile(d(idx_windows).folder, d(idx_windows).name, sprintf('fit_results_%s.mat', model_str)));
    plot_estimates('base', fh, 'results', model_results, 'marker', markers{i}, 'label', d(idx_windows).name)
    legend('Location', 'eastoutside', 'box', 'off', 'FontSize', 20)

end