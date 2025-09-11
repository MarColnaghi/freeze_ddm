function [lls_table] = test_consistency(run_code, model, varargin)

opt = inputParser;
addParameter(opt, 'export', false);
addParameter(opt, 'type', 'surrogate');
addParameter(opt, 'gen_model', '');  % Plot all trials if true

parse(opt, varargin{:});

export = opt.Results.export;
type = opt.Results.type;
gen_model = opt.Results.gen_model;

col = cmapper();
lls_table = table();

id = 0;
ll_id = table();
ll_tv = table();
ll_cf = table();
ll_st = table();
lls_table = table();

for idx_run = run_code
    
    id = id + 1;

    run_str = sprintf('run%02d', idx_run);
    run_path = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/momentary_integration', type, model, run_str);
    cd(run_path);

    load('rt.mat')
    list = dir(fullfile(run_path, sprintf('*%s/lls_final.mat', gen_model)));
    
    for idx_file = 1:length(list)
        fullfile(list(idx_file).folder, list(idx_file).name)
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

    end

    ll_id.(run_str) = lls.id; 
    ll_tv.(run_str) = lls.ll_tv; 
    ll_st.(run_str) = lls.ll_st; 
    ll_cf.(run_str) = lls.ll_cf; 

end

for idx_gen_model = {'id', 'tv', 'st', 'cf'}
    lls_table.(idx_gen_model{1}) = eval(sprintf('ll_%s', idx_gen_model{1}));
end

fh = figure('color','w','Position',[100, 100, 950, 400]);
tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
hold on;
i = 0;
for idx_gen_model = {'tv', 'st', 'cf'}
    nexttile
    i = i + 1;
    disp(idx_run)
    corrx = corr(table2array(lls_table.(idx_gen_model{1})));
    
    imagesc(corrx, [0 1])

    apply_generic(gca)
    axis square
    xticks(1:size(corrx,2))
    yticks(1:size(corrx,2))
    
end
colorbar(gca, 'eastoutside', 'LineWidth', 2, 'FontSize', 20)
if export
    cd(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/figures/momentary_integration', type, model));
    paths.fig = fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/figures/momentary_integration', type);
    exporter(fh, paths, 'consistency.pdf');
end
end