function analyse_evidence_traj(run_code, model, varargin)

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
        all_traj = importdata(fullfile(list(idx_file).folder, 'all_traj.mat'));
        cell_sm = extract_sm_cellarray(y, motion_cache);

        [log, idx] = ismember(lls.id, rt.id);
        sort_rt = rt(idx,:);

        if any(~log)
            disp('missing some datapoints')
        end
        if sort_rt.id(1) == lls.id(1)
            disp('matching')
        end

        delta_ll = lls.ll_tv - lls.ll_st;
        [sorted_deltall, idx_deltall] = sort(delta_ll);
        lls
        for idx_bouts = 1:size(all_traj)
            sm = cell_sm{idx_bouts};
            evidence = all_traj{idx_bouts};
            evidence = evidence(~isnan(evidence));
            outcorr = corrcoef(evidence, sm(1:length(evidence)));
            corr(idx_bouts) = outcorr(2);

%             if abs(corr(idx_bouts)) > 0.7 && y.durations_s(idx_bouts) > 0.5
%                 figure
%                 plot(evidence)
%                 hold on
%                 plot(sm(1:length(evidence)))
%                 text(1,1,num2str(corr(idx_bouts)))
%                 text(1,0.8,num2str(delta_ll(idx_bouts)))
%                 pause(0.1)
%             end
        end

        figure
        histogram(corr, [-1:0.025:1])
        figure
        hold on
        mask = y.durations_s > 5;
        scatter(corr, delta_ll, 5, 'k')
        scatter(corr(mask), delta_ll(mask), 5, 'r')
        fh = figure('color','w','Position',[100, 100, 1600, 800]);
        n_col = 3; n_rows = 3; items = n_col * n_rows;
        tiledlayout(n_rows, n_col, 'TileSpacing', 'tight', 'Padding', 'compact')
        
        i = 0;
        if strcmp(selection, 'best')
            selected_trials = idx_deltall(end - items + 1:end);
        elseif strcmp(selection, 'worst')
            selected_trials = idx_deltall(1:items);
        end

        for idx_tile = 1:items
            nexttile
            plot(all_traj{selected_trials(idx_tile)});

            hold on
            plot(cell_sm{selected_trials(idx_tile)});
            xlim([0 660])
            ylim([-5 10])
            apply_generic(gca)
        end
        if export
            cd(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/figures/momentary_integration', type, model, run_str));
            folderi = sprintf('sims_%s', gen_data);
            mkdir(folderi); paths.fig = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/figures/momentary_integration', type, model, run_str, folderi);
            exporter(fh, paths, sprintf('examples_%s.pdf', selection));
        end
    end
end

function apply_legend(ax)

ax(1).XAxis.TickValues = [0 15];
xlabel('')
