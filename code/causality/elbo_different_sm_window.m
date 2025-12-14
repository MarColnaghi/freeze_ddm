% Double DDM test

clearvars

% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 2; threshold_mob = 2; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
model_str = 'dddm2';
col.Set3 = cbrewer2('Set3', 5);

windows.anchor = 'loom_onset';
windows.reference = 'fixed_length';
windows.length = '30';

fh = figure('Position', [100 100 1070 500], 'Color', 'w');
t = tiledlayout(1, 1, 'TileSpacing', 'loose', 'Padding', 'loose');
mkrsize = 180;

model_acausal = importdata('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/fitting_freezes/le/dddm2/run11/fit_results_dddm2.mat');
elbo_acausal = model_acausal.elbo;

parent_folder = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/causality/fitting_windows', model_str, windows.reference);
code_folder = sprintf('*_size%s*', windows.length);
dir_folder = dir(fullfile(parent_folder, code_folder));

nexttile

c = 0;
for idx_run = 1:length(dir_folder)
    
    c = c + 1;
    str_folder = fullfile(dir_folder(idx_run).folder, dir_folder(idx_run).name);

    extractLabel = @(s) ...
        regexprep( ...
        regexprep(extractAfter(s, '-'), '_', ' '), ...
        '(^|\s)(\w)', '${upper($2)}');

    label = extractLabel(dir_folder(idx_run).name);

    % paths = path_generator('folder', 'causality/fixed_window', 'bouts_id', id_code, 'imfirst', false);
    d = dir(str_folder);
    d = d(~ismember({d.name}, {'.','..', '.DS_Store'}));

    elbos = nan(size(d, 1), 2);
    get_window_size = nan(size(d, 1), 1);

    for idx_windows = 1:size(d, 1)
        get_window_size(idx_windows) = str2double(d(idx_windows).name);
        model_results = importdata(fullfile(d(idx_windows).folder, d(idx_windows).name, sprintf('fit_results_%s.mat', model_str)));
        elbos(idx_windows, :) = model_results.elbo;
        
    end

    hold on

    scatter(get_window_size, elbos(:,1), mkrsize, 'o', 'MarkerFaceColor', col.Set3(3 + c, :), ...
        'MarkerEdgeColor', 'none', 'DisplayName', sprintf('Aligned to %s', label))
    [sorted_window_size, i] = sort(get_window_size);
    plot(sorted_window_size, elbos(i, 1), 'Color', col.Set3(3 + c, :), 'LineStyle', '--', 'HandleVisibility','off')
    grid on
end



% scatter(30, elbo_acausal(:,1), mkrsize, 'o', 'MarkerFaceColor', 'k', ...
%         'MarkerEdgeColor', 'none')

yline(elbo_acausal(:,1), 'k--', 'DisplayName', 'Freeze Duration', 'LineWidth', 2)
ax = gca;
apply_generic(ax, 'xlim', [-610 240], 'xticks', [-600 -300 -60 0 30 60 120 180 240], 'ylim', [-10650 -10400])
ax.GridAlpha = .1;

if strcmp(windows.reference, 'fixed_length')
    xlabel('Start of Integration Windows (frames)')
else
    xlabel('Start/End of Integration Window (frames)')
end

ylabel('ELBO')
legend('Location', 'southeastoutside', 'Box', 'off', 'FontSize', 20)

if strcmp(windows.reference, 'fixed_length')
    figure_label = sprintf('elbos_%s_length%s.pdf', windows.reference, windows.length);
else
    figure_label = sprintf('elbos_%s.pdf', windows.reference);
end

paths = path_generator('folder', fullfile('causality/fitting_windows', model_str, windows.reference));
exporter(fh, paths, figure_label)