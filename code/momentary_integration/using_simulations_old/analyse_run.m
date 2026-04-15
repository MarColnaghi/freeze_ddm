function analyse_run(run_code, model, varargin)

opt = inputParser;
addParameter(opt, 'export', false);
addParameter(opt, 'type', 'surrogate');
addParameter(opt, 'gen_model', '');  % Plot all trials if true

parse(opt, varargin{:});

export = opt.Results.export;
type = opt.Results.type;
gen_model = opt.Results.gen_model;

col = cmapper();


for idx_run = run_code

    run_str = sprintf('run%02d', idx_run);
    cd(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/momentary_integration', type, model, run_str))
    
    load('rt.mat')
    list = dir(sprintf('*%s/lls_final.mat', gen_model));
    
    for idx_file = 1:length(list)
        lls = importdata(fullfile(list(idx_file).folder, list(idx_file).name));
        gen_data = list(idx_file).folder(end-1:end);

        [log, idx] = ismember(lls.id, rt.id);
        sort_rt = rt(idx,:);

        if any(~log)
            disp('missing some datapoints')
        end
        if sort_rt.id(1) == lls.id(1)
            disp('matching')
        end

        fh = figure('color','w','Position',[100, 100, 800, 300]);
        tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact')
        nexttile
        total_ll_bar(lls, col, gen_data)
        
        nexttile(2, [1 2])
        hold on
        ylim([-3 3])
        xlim([-50 height(lls) + 50])
        ax = gca;

        [sorted_deltall, idx_deltall_tv] = sort(lls.ll_tv - lls.ll_st);
        idx_crossing = find(sorted_deltall > 0);

        fill([ax.XLim(1) ax.XLim(1) idx_crossing(1) idx_crossing(1)], [ax.YLim(1) 0 0 ax.YLim(1)], hex2rgb(col.stationary_sm), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'LineWidth', 2)
        fill([idx_crossing(1) idx_crossing(1) ax.XLim(2) ax.XLim(2)], [ax.YLim(2) 0 0 ax.YLim(2)], hex2rgb(col.timevarying_sm), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'LineWidth', 2)

        xline(idx_crossing(1), 'Label', sprintf('%.2f%%', (sum(sorted_deltall > 0)./length(sorted_deltall)) * 100), 'FontSize', 20, 'LabelOrientation','horizontal')

        bar(sorted_deltall, 'FaceColor', 'k', 'EdgeColor', 'none' )
        xlabel('Sorted Freezes')
        xticks([])
        ylabel('$\log \!\big(\mathcal{L}_{\mathrm{tv}} / \mathcal{L}_{\mathrm{st}}\big)$', ...
       'Interpreter','latex')
        apply_generic(ax)

        if export
            cd(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/figures/momentary_integration', type, model, run_str));
            folderi = sprintf('sims_%s', gen_data);
            mkdir(folderi); paths.fig = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/figures/momentary_integration', type, model, run_str, folderi);
            exporter(fh, paths, 'total_ll.eps')
        end
    end
end
end

function total_ll_bar(lls_output, col, gen_data)

bh = bar([sum(lls_output.ll_tv), sum(lls_output.ll_st), sum(lls_output.ll_cf)], 'FaceColor', 'flat', 'EdgeColor', 'flat', 'LineWidth', 2);
bh.CData(1,:) = hex2rgb(col.timevarying_sm);
bh.CData(2,:) = hex2rgb(col.stationary_sm);
bh.CData(3,:) = hex2rgb(col.stationary_sm);

if strcmp(gen_data, 'tv')
        text(1, bh.YData(1), 'generative', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 20, 'Color', 'w', 'Rotation', 90)

elseif strcmp(gen_data, 'st')
        text(2, bh.YData(2), 'generative', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 20, 'Color', 'w', 'Rotation', 90)

elseif strcmp(gen_data, 'ig')
        text(3, bh.YData(3), 'generative', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', 'FontSize', 20, 'Color', 'w', 'Rotation', 90)

end

xticklabels({'tv', 'st', 'st-cf'})
ylabel('Total log($\mathcal{L}$)', 'Interpreter', 'latex')
ax = gca;
apply_generic(ax)
ax.YAxis.Direction = 'reverse';
ylim([-15000 0]); %ax.YTickLabel(1) = {'worse'}; ax.YTickLabel(end) = {'better'};
text(2, ax.YLim(1) + 1000, sprintf('$\\Delta log(\\mathcal{L}): %.2f$', sum(lls_output.ll_tv) - sum(lls_output.ll_st)), 'FontSize', 20, 'Interpreter', 'latex', 'HorizontalAlignment', 'center')
end