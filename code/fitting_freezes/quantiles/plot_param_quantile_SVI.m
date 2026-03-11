paths = path_generator('folder', 'fitting_freezes/le/quantiles', 'bouts_id', id_code);

imported_results = importdata('/Users/marcocolnaghi/PhD/bayes_fpe/model_results_WO_FS/Master_Summary_SVI.csv');
indices_mean = find(contains(imported_results.textdata, 'mean', 'IgnoreCase', true));
indices_std = find(contains(imported_results.textdata, 'std', 'IgnoreCase', true));

% 1. Define the exact order you want
meanLabels = {'drift_1_mean', 'bound_1_mean', 'drift_2_mean', 'bound_2_mean', 'pmix_mean', 'shared_ndt_mean'};
stdLabels = {'drift_1_std', 'bound_1_std', 'drift_2_std', 'bound_2_std', 'pmix_std', 'shared_ndt_std'};
% Autocompleted labels with LaTeX math symbols
titleLabels = {...
    '$\mathrm{Drift}_{1} (\mu_{1})$', ...
    '$\mathrm{Bound}_{1} (\theta_{1})$', ...
    '$\mathrm{Drift}_{2} (\mu_{2})$', ...
    '$\mathrm{Bound}_{2} (\theta_{2})$', ...
    '$p_{mix}$', ...
    '$ndt$'};


% 2. Find the column indices (sourceIdx)
[~, means] = ismember(targetLabels, imported_results.textdata);
[~, stds] = ismember(stdLabels, imported_results.textdata);

% 3. Extract the numeric columns from the data matrix
% This takes all rows (:) and only the columns matching our targets
meanValues = imported_results.data(:, means);
stdValues = imported_results.data(:, stds);

% 4. (Optional) Create a Table for easy viewing
meanTable = array2table(meanValues, 'VariableNames', {'drift_1', 'bound_1', 'drift_2', 'bound_2', 'pmix', 'shared_ndt'});
stdTable = array2table(stdValues, 'VariableNames', {'drift_1', 'bound_1', 'drift_2', 'bound_2', 'pmix', 'shared_ndt'});


fh = figure('color','w','Position',[100,100, 1200, 500]);
tiledlayout(1, 6, 'TileSpacing', 'loose', 'Padding', 'compact')

for idx_params = 1:6
    nexttile
    hold on
%     scatter([1:3] - 0.05, table2array(meanTable(imported_results.data(:,3) == 25, idx_params)), 50, col.vars.ls(3,:), 'filled')
%     scatter([1:3] + 0.05, table2array(meanTable(imported_results.data(:,3) == 50, idx_params)), 50, col.vars.ls(5,:), 'filled')
%     errbar([1:3] - 0.05,  table2array(meanTable(imported_results.data(:,3) == 25, idx_params)), ....
%         table2array(stdTable(imported_results.data(:,3) == 25, idx_params)),'Color', col.vars.ls(3,:))
%     errbar([1:3] + 0.05,  table2array(meanTable(imported_results.data(:,3) == 50, idx_params)), ....
%         table2array(stdTable(imported_results.data(:,3) == 50, idx_params)), 'Color', col.vars.ls(5,:))
     axis square


    scatter([1:3] - 0.05, table2array(meanTable(imported_results.data(:,3) == 25, idx_params)), 50, [0.4 0.4 0.4], 'filled')
%     scatter([1:3] + 0.05, table2array(meanTable(imported_results.data(:,3) == 50, idx_params)), 50, col.vars.ls(5,:), 'filled')
    errbar([1:3] - 0.05,  table2array(meanTable(imported_results.data(:,3) == 25, idx_params)), ....
        table2array(stdTable(imported_results.data(:,3) == 25, idx_params)),'Color', [0.4 0.4 0.4])
%     errbar([1:3] + 0.05,  table2array(meanTable(imported_results.data(:,3) == 50, idx_params)), ....
%         table2array(stdTable(imported_results.data(:,3) == 50, idx_params)), 'Color', col.vars.ls(5,:))

    xlim([0.5 3.5])
    xticks(1:3)
    xticklabels([1 2 3 4])

    if idx_params == 1
        ylim([0 3])
    elseif idx_params == 2
        ylim([0 1])

    elseif idx_params == 3
        ylim([0 1.5])

    elseif idx_params == 4
        ylim([0 5])

    elseif idx_params == 5
        ylim([0 1])

    elseif idx_params == 6
        ylim([0 .5])
    end

    apply_generic(gca, 'line_width', 2, 'tick_length', 0.02)

    xlabel({'Soc. Mot.', 'Range'}, 'FontSize', 18)
    ax = gca;
    ax.XAxis.FontSize = 18;
    xticklabels({'Low', 'Med', 'High'})

    title(titleLabels{idx_params}, 'Interpreter', 'latex', 'FontSize', 22);
end

exporter(fh, paths, 'quartile_fitting_SVI_onlyslow.pdf')
