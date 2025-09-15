function [early_sm, late_sm, com, slopes, Signal_Slope, Center_of_Mass, Trend_Score, ES_difference, Variance, sort_ll] = analyse_sm_detailed(run_code, model, varargin)
    % ANALYSE_SM Visualize motion traces sorted by model delta-LLs
    % Usage: analyse_sm(run_code, model, 'type', 'surrogate', 'selection', 500, 'export', true, 'all', false)

    %% Parse inputs
    p = inputParser;
    addParameter(p, 'export', false);
    addParameter(p, 'type', 'surrogate');
    addParameter(p, 'selection', 500);
    addParameter(p, 'gen_model', '');  % Plot all trials if true
    addParameter(p, 'all', false);  % Plot all trials if true
    addParameter(p, 'align_to', 'offset');  % 'offset' (old behavior) or 'onset'

   
    parse(p, varargin{:});

    export = p.Results.export;
    type = p.Results.type;
    selection = p.Results.selection;
    plot_all = p.Results.all;
    gen_model = p.Results.gen_model;
    align_to = p.Results.align_to;

    col = cmapper();  % Your custom colormap function
    close all
    
    min_dur = 0.35;
    min_sm = 0.1;
    max_sm = 0.7;
    opt.thr_slope = 0.0;
    opt.thr_com = 0.42; % Set 1 if you don't want it
    opt.mov_mean = 1; % Set 1 if you don't want it
    opt.thr_trend = 0.55; % Set 1 if you don't want it
    opt.thr_es = .8;

    for idx_run = run_code
        run_str = sprintf('run%02d', idx_run);
        run_path = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/momentary_integration', type, model, run_str);
        fprintf('Processing %s...\n', run_path);

        %% Load necessary data
        load(fullfile(run_path, 'model_results.mat'), 'model_results');
        load(fullfile(run_path, 'sim_params.mat'), 'sim_params');
        load(fullfile(run_path, 'rt.mat'), 'rt');
        load(sim_params.motion_cache_path, 'motion_cache');
        colmap = cbrewer2('RdbU', height(rt));
        colmap = horzcat(colmap, 0.1 * ones(height(rt),1));

        list = dir(fullfile(run_path, sprintf('*%s/lls_final.mat', gen_model)));

        for idx_file = 1:length(list)
            folder = list(idx_file).folder;
            gen_data = folder(end-1:end);
            
            lls = importdata(fullfile(folder, 'lls_final.mat'));
            load(fullfile(folder, 'y.mat'));

            [~, idx_rt] = ismember(lls.id, rt.id);
            if any(idx_rt == 0)
                warning('Some lls IDs not found in RT.');
                continue;
            end
            
            sort_rt = rt(idx_rt, :);

            %% Sanity check
            if isequal(sort_rt.id(1), lls.id(1))
                disp('Matching IDs confirmed');
            end

            %% Sort by delta LL
            delta_ll = lls.ll_tv - lls.ll_st;
            [sort_ll, sort_idx] = sort(delta_ll);
            y_sorted = y(sort_idx, :);

            %% Process traces
            trace_len = length(motion_cache(1));
            cell_sm = extract_sm_cellarray(y_sorted, motion_cache);
            [early_sm, late_sm, com, slopes, Signal_Slope, Center_of_Mass, Trend_Score, ES_difference, Variance] =  partition_variable_freezes(cell_sm, opt);



            
            metrics = {'Signal_Slope', 'Center_of_Mass', 'ES_difference', 'Variance'};
            limits = {[-0.1 0.1], [0 1], [-2 2], [0 1], [0.45 0.55]};

            for idx_metrics = 1:length(metrics)
                fh = figure('color', 'w', 'Position', [100, 100, 800, 350]);
                tl = tiledlayout(1, 3, 'TileSpacing', 'loose', 'Padding', 'loose');
                nexttile
                %imagesc(eval(metrics{idx_metrics}), limits{idx_metrics})
                scatter(1:height(y_sorted), eval(metrics{idx_metrics}), 5, 'k', 'filled')
                apply_generic(gca, 20)
                xlim([-50 height(y_sorted) + 50])
                ylim(limits{idx_metrics})
                ax = gca;
                xticks([]); xlabel('sorted'); ylabel(strrep(metrics{idx_metrics}, '_', ' '));

                nexttile
                scatter(sort_ll, eval(metrics{idx_metrics}), 5, 'k', 'filled')
                ylim(limits{idx_metrics})
                ax = gca;
                apply_generic(ax)
                ax.YAxis.Visible = 'off';
                xlim([-3 3]); xlabel('$\Delta_{ll}$', 'Interpreter', 'latex');

                nexttile
                scatter(sort_ll, eval(metrics{idx_metrics}), 5, y_sorted.durations_s, 'filled')
                clim([0 5])
                colorcet('I1')
                ax = gca;
                apply_generic(ax)
                ylim(limits{idx_metrics})
                xlim([-3 3]); %xlabel('$\Delta_{ll}$', 'Interpreter', 'latex');
                ax.YAxis.Visible = 'off';
                cb = colorbar(ax, 'location', "eastoutside", "LineWidth", 2, 'FontSize', 20);
                cb.Label.String = 'Freeze Duration (s)';
                xline(0)
                yline(sum(limits{idx_metrics})./2)
                
                if export
                    cd(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/figures/momentary_integration', type, model, run_str));
                    folderi = sprintf('sims_%s', gen_data);
                    mkdir(folderi); paths.fig = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/figures/momentary_integration', type, model, run_str, folderi);
                    exporter(fh, paths, sprintf('analyse_sm_%s.pdf', metrics{idx_metrics}));
                end


            end


        end
    end
end
