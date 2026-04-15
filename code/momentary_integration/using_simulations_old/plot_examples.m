function plot_examples(run_code, model, varargin)

opt = inputParser;
addParameter(opt, 'export', false);
addParameter(opt, 'type', 'surrogate');
addParameter(opt, 'selection', 'best');
addParameter(opt, 'gen_model', '');  % Plot all trials if true

parse(opt, varargin{:});

export = opt.Results.export;
type = opt.Results.type;
selection = opt.Results.selection;
gen_model = opt.Results.gen_model;
col = cmapper();
sim_params.num_sims = 100000;

for idx_run = run_code

    run_str = sprintf('run%02d', idx_run);
    run_path = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/momentary_integration', type, model, run_str);
    cd(run_path)


    load('model_results.mat');
    load('sim_params.mat');
    load('rt.mat'); 
    load(sim_params.motion_cache_path)

    list = dir(fullfile(run_path, sprintf('*%s/lls_final.mat', gen_model)));

    for idx_file = 1:length(list)
        lls = importdata(fullfile(list(idx_file).folder, list(idx_file).name));
        gen_data = list(idx_file).folder(end-1:end);
        load(fullfile(list(idx_file).folder, 'y.mat'));

        [log, idx] = ismember(lls.id, rt.id);
        sort_rt = rt(idx,:);

        if any(~log)
            disp('missing some datapoints')
        end
        if sort_rt.id(1) == lls.id(1)
            disp('matching')
        end

        [deltall, idx_deltall] = sort(lls.ll_tv - lls.ll_st);
    
        fh = figure('color','w','Position',[100, 100, 1200, 800]);
        n_col = 3; n_rows = 3; items = n_col * n_rows;
        tiledlayout(n_rows, n_col, 'TileSpacing', 'tight', 'Padding', 'compact')
        
        i = 0;
        if strcmp(selection, 'best')
            selected_trials = idx_deltall(end - items + 1:end);
        elseif strcmp(selection, 'worst')
            selected_trials = idx_deltall(1:items);
        elseif strcmp(selection, 'any')
            selected_trials = randi(length(idx_deltall), items, 1);

        end

        for idx_tile = 1:items
            nexttile
            if idx_tile == 1
                evaluate_likelihood_plot(y, motion_cache, sim_params, model_results.(gen_data), selected_trials(idx_tile), 'legend', true);
            else
                evaluate_likelihood_plot(y, motion_cache, sim_params, model_results.(gen_data), selected_trials(idx_tile), 'legend', false);

            end

        end
        if export
            cd(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/figures/momentary_integration', type, model, run_str));
            folderi = sprintf('sims_%s', gen_data);
            mkdir(folderi); paths.fig = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/figures/momentary_integration', type, model, run_str, folderi);
            exporter(fh, paths, sprintf('examples_%s.eps', selection));
        end
    end
end


function apply_legend(ax)

ax(1).XAxis.TickValues = [0 11];
xticks([0 10])
xlabel('')
