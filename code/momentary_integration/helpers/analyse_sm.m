function analyse_sm(run_code, model, varargin)
    % ANALYSE_SM Visualize motion traces sorted by model delta-LLs
    % Usage: analyse_sm(run_code, model, 'type', 'surrogate', 'selection', 500, 'export', true, 'all', false)

    %% Parse inputs
    p = inputParser;
    addParameter(p, 'export', false);
    addParameter(p, 'type', 'surrogate');
    addParameter(p, 'selection', 500);
    addParameter(p, 'gen_model', '');  % Plot all trials if true
    addParameter(p, 'all', false);  % Plot all trials if true

   
    parse(p, varargin{:});

    export = p.Results.export;
    type = p.Results.type;
    selection = p.Results.selection;
    plot_all = p.Results.all;
    gen_model = p.Results.gen_model;

    col = cmapper();  % Your custom colormap function
    
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
            gen_data = folder(end-1:end)
            
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
            [~, sort_idx] = sort(delta_ll);
            y_sorted = y(sort_idx, :);

            %timelock_sm
            %[sorted_len, sort_idx] = sort(cellfun(@length, timelock_sm), 'descend');
            %aligned_mat = padder_for_imagesc(timelock_sm(sort_idx), sorted_len, 'offset');

            %% Process traces
            trace_len = length(motion_cache(1));

            if plot_all
                % Store all trials
                all_mat = nan(height(y_sorted), trace_len);
                for i = 1:height(y_sorted)
                    trace = get_cropped_signal(motion_cache, y_sorted, i, trace_len, model_results.(gen_data).estimates_mean);
                    all_mat(i, :) = trace;
                end
            else
                % Top and bottom trial selection
                num_trials = selection;
                top_mat = nan(num_trials, trace_len);
                bottom_mat = nan(num_trials, trace_len);

                for i = 1:num_trials
                    trace = get_cropped_signal(motion_cache, y_sorted, i, trace_len, model_results.(gen_data).estimates_mean);
                    bottom_mat(i, :) = trace;
                end

                for i = height(y_sorted) - num_trials + 1 : height(y_sorted)
                    trace = get_cropped_signal(motion_cache, y_sorted, i, trace_len, model_results.(gen_data).estimates_mean);
                    top_mat(i - height(y_sorted) + num_trials, :) = trace;
                end
            end

            %% Plot
            fh = figure('color','w','Position',[100, 100, 600, 400]);
            hold on;

            if plot_all
                for idx_trials = 1:height(lls)
                    plot(- trace_len + 1:1:0, all_mat(idx_trials,:), 'Color', colmap(idx_trials,:),'LineWidth', 0.25);
                end
            else
                plot(-trace_len + 1:1:0, bottom_mat, 'Color', [hex2rgb(col.stationary_sm) 0.05], 'LineWidth', 0.5);
                plot(-trace_len + 1:1:0, median(bottom_mat, 'omitnan'), 'Color', col.stationary_sm, 'LineWidth', 2);
                plot(-trace_len + 1:1:0, top_mat, 'Color', [hex2rgb(col.timevarying_sm) 0.05], 'LineWidth', 0.5);
                plot(-trace_len + 1:1:0, median(top_mat, 'omitnan'), 'Color', col.timevarying_sm, 'LineWidth', 2);
            end

            xlabel('Time to unfreeze (frames)');
            ylabel('Motion Signal (relative)');
            apply_generic(gca);

            xlim([-600 0]);
            ylim([-3 3]);

            if export
                cd(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/figures/momentary_integration', type, model, run_str));
                folderi = sprintf('sims_%s', gen_data);
                mkdir(folderi); paths.fig = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/figures/momentary_integration', type, model, run_str, folderi);
                exporter(fh, paths, 'analyse_sm.pdf');
            end
        end
    end
end

%% Helper function
function cropped = get_cropped_signal(motion_cache, y, idx, trace_len, estimates_mean)
    signal = motion_cache(y.fly(idx));
    cropped = nan(1, trace_len);
    onset = y.onsets(idx);
    duration_frames = round(y.durations_s(idx) * 60 - estimates_mean.tndt_intercept * 60);

    end_idx = onset + duration_frames - 1;
    if onset < 1 || end_idx > length(signal)
        return;
    end

    cropped(end - duration_frames + 1:end) = signal(onset:end_idx);
    cropped = cropped - cropped(end);  % Align endpoint to 0
end
