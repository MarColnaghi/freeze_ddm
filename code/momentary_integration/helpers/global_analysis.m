function global_analysis(run_code, model, varargin)

opt = inputParser;
addParameter(opt, 'export', false);
addParameter(opt, 'type', 'surrogate');
addParameter(opt, 'gen_model', '');  % Plot all trials if true

parse(opt, varargin{:});

export = opt.Results.export;
type = opt.Results.type;
gen_model = opt.Results.gen_model;

col = cmapper();


fh = figure('color','w','Position',[100, 100, 400, 200]);
hold on

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
        
        if sum(lls.ll_tv) - sum(lls.ll_st) > 0
            scatter(1, sum(lls.ll_tv) - sum(lls.ll_st), 120, hex2rgb(col.timevarying_sm), 'filled', 'MarkerFaceAlpha', 0.35)
        elseif sum(lls.ll_tv) - sum(lls.ll_st) < 0
            scatter(2, sum(lls.ll_tv) - sum(lls.ll_st), 120, hex2rgb(col.stationary_sm), 'filled', 'MarkerFaceAlpha', 0.35)
        end

    end
end

for idx_run = 1:5
    run_str = sprintf('run%02d', idx_run);
    cd(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/momentary_integration', 'freezes', model, run_str))
    
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

        scatter(1.5, sum(lls.ll_tv) - sum(lls.ll_st), 120, 'k', 'filled', 'MarkerFaceAlpha', 0.35)


    end
end


apply_generic(gca)
xlim([0.5 2.5])
xticks([1,1.5,2])
xticklabels({'sim. tv', 'data', 'sim. st'})
xtickangle(0)
ylim([-500 500])
%ylabel('delta LL')
ax = gca;

if export
    cd(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/figures/momentary_integration', 'freezes', model));
    paths.fig = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/figures/momentary_integration', 'freezes', model);
    exporter(fh, paths, 'summary.pdf')
end

end
