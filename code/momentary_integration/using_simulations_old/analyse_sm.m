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
    addParameter(p, 'align_to', 'offset');  % 'offset' (old behavior) or 'onset'

   
    parse(p, varargin{:});

    export = p.Results.export;
    type = p.Results.type;
    selection = p.Results.selection;
    plot_all = p.Results.all;
    gen_model = p.Results.gen_model;
    align_to = p.Results.align_to;

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
                for i = 1:height(y_sorted)
                    trace = get_cropped_signal(motion_cache, y_sorted, i, trace_len, model_results.(gen_data).estimates_mean, align_to);
                    all_mat(i, :) = trace;
                end
            else
                for i = 1:selection
                    trace = get_cropped_signal(motion_cache, y_sorted, i, trace_len, model_results.(gen_data).estimates_mean, align_to);
                    bottom_mat(i, :) = trace;
                end
                for i = height(y_sorted) - selection + 1 : height(y_sorted)
                    trace = get_cropped_signal(motion_cache, y_sorted, i, trace_len, model_results.(gen_data).estimates_mean, align_to);
                    top_mat(i - height(y_sorted) + selection, :) = trace;
                end
            end

            % ----- X axis depends on alignment -----
            if strcmpi(align_to,'onset')
                x = 0:trace_len-1;                  % frames since onset
                xlabel({'Time (frames)', '[freeze start aligned]'});
                ylabel({'Soc. Motion Signal', '(start subtracted)'});
            else
                x = -trace_len+1:0;                  % old behavior: time to unfreeze
                xlabel({'Time (frames)', '[freeze break aligned]'});
                ylabel({'Soc. Motion Signal', '(endpoint subtracted)'});
            end

            % ----- Plot -----
            fh = figure('color','w','Position',[100, 100, 600, 400]); hold on;
            if plot_all
                for idx_trials = 1:height(lls)
                    plot(x, all_mat(idx_trials,:), 'Color', colmap(idx_trials,:),'LineWidth', 0.25);
                end
            else
                plot(x, bottom_mat, 'Color', [hex2rgb(col.stationary_sm) 0.1], 'LineWidth', 0.5);
                plot(x, median(bottom_mat, 1, 'omitnan'), 'Color', col.stationary_sm, 'LineWidth', 2);
                plot(x, top_mat, 'Color', [hex2rgb(col.timevarying_sm) 0.1], 'LineWidth', 0.5);
                plot(x, median(top_mat, 1, 'omitnan'), 'Color', col.timevarying_sm, 'LineWidth', 2);
            end

            apply_generic(gca);

            if strcmpi(align_to, 'onset')
                xlim([0 600]);   % adjust as you like
                xlabel({'Time (frames)', '[freeze start aligned]'});
                ylabel({'Soc. Motion Signal', '(start subtracted)'});

            else
                xlim([-600 0]);
                xlabel({'Time (frames)', '[freeze break aligned]'});
                ylabel({'Soc. Motion Signal', '(endpoint subtracted)'});

            end
            ylim([-3 3]);

            if export
                cd(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/figures/momentary_integration', type, model, run_str));
                folderi = sprintf('sims_%s', gen_data);
                mkdir(folderi); paths.fig = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/figures/momentary_integration', type, model, run_str, folderi);
                exporter(fh, paths, sprintf('analyse_sm_%s_%d.pdf', align_to, selection));
            end
        end
    end
end

%% Helper function
function cropped = get_cropped_signal(motion_cache, y, idx, trace_len, estimates_mean, align_to)
    if nargin < 6 || isempty(align_to); align_to = 'offset'; end

    signal = motion_cache(y.fly(idx));
    cropped = nan(1, trace_len);

    onset = y.onsets(idx);
    duration_frames = round(y.durations_s(idx) * 60 - estimates_mean.tndt_intercept * 60);
    if duration_frames < 1; return; end

    end_idx = onset + duration_frames - 1;
    if onset < 1 || end_idx > numel(signal)
        return; % out-of-bounds, keep NaNs
    end

    seg = signal(onset:end_idx);

    if strcmpi(align_to,'onset')
        % Place segment at the BEGINNING of the buffer
        n = min(numel(seg), trace_len);
        cropped(1:n) = seg(1:n);
        % Zero at onset (first valid sample)
        if ~isnan(cropped(1)); cropped = cropped - cropped(1); end
    else
        % Old behavior: place segment at the END of the buffer
        n = min(numel(seg), trace_len);
        cropped(end-n+1:end) = seg(end-n+1:end);
        % Zero at the endpoint (last valid sample)
        if ~isnan(cropped(end)); cropped = cropped - cropped(end); end
    end
end

