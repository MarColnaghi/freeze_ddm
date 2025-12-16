function [fh] = base_for_estimates(varargin)

opt = inputParser;

addParameter(opt, 'model', []);          % <-- NEW
addParameter(opt, 'export', false);
addParameter(opt, 'ylimits', [-1, 4]);
addParameter(opt, 'paths', []);
addParameter(opt, 'fh', []);

% optional: pass estimates if you want later
addParameter(opt, 'est_means', []);      % numeric vector or table

% optional: keep order stable (because fieldnames() is alphabetical)
addParameter(opt, 'component_order', {});  % e.g. {'mu1','theta1','mu2','theta2','pmix','tndt'}

parse(opt, varargin{:});
model = opt.Results.model;
export = opt.Results.export; %#ok<NASGU>
paths  = opt.Results.paths;  %#ok<NASGU>
ylimits = opt.Results.ylimits;
fh = opt.Results.fh;

% ---- Flatten model -> parameter table with VariableNames like mu1_sm, theta2_intercept, ...
param_tbl = flatten_model_to_param_table(model, opt.Results.component_order);

% If you want to keep compatibility with your old flow:
% "est_means" here is just a table with variable names; values are optional.
est_means = param_tbl;

quantiles = 0;
col = cmapper([], quantiles);
if isempty(fh)
    fh = figure('color','w', 'Position', [100 100 800 400]);
end
hold on

[suffixes, prefixes] = extract_dep(est_means);

% Replace 'intercept' with '0' in the suffixes (your existing logic)
suffixes_replaced = suffixes;
suffixes_replaced(strcmp(suffixes, 'intercept')) = {'0'};

xx = 1:numel(suffixes);

% If you want: β_{component}^{predictor}
result = cellfun(@(pre, suf) ['$$\beta_{' pre '}^{' suf '}$$'], ...
                 prefixes, suffixes_replaced, 'UniformOutput', false);

% (You had a second overwrite of `result` in your snippet; removing that to keep the intended label.)

for idx_param = 1:length(xx)
    % background patch per parameter
    this_suf = suffixes{idx_param};

    % guard: if cmapper doesn't have this field, fallback to gray-ish
    if isfield(col.vars, this_suf)
        faceCol = col.vars.(this_suf);
    else
        faceCol = [0.5 0.5 0.5];
    end

    fill([xx(idx_param) - 0.3, xx(idx_param) - 0.3, xx(idx_param) + 0.3, xx(idx_param) + 0.3], ...
         [ylimits, fliplr(ylimits)], '', ...
         'FaceColor', faceCol, 'LineStyle', 'none', 'FaceAlpha', 0.14, 'HandleVisibility', 'off');
end

xlim([xx(1) - 1, xx(end) + 1]);
ylim(ylimits);

ax = gca;
apply_generic(ax)

xticks(xx);
xticklabels(result);
set(ax.XAxis, 'TickLabelInterpreter', 'latex', 'FontSize', 24);

end


% ---------- helper: flatten model struct to a "fake estimates table" ----------
function param_tbl = flatten_model_to_param_table(model, component_order)

if isempty(model) || ~isstruct(model)
    error('base_for_estimates:InvalidModel', 'You must pass a model struct (e.g. model_dddm2()).');
end

% Choose component order
if ~isempty(component_order)
    comps = component_order;
else
    comps = fieldnames(model); % NOTE: alphabetical by default
end

varNames = {};

for i = 1:numel(comps)
    comp = comps{i};
    if ~isfield(model, comp), continue; end

    block = model.(comp);
    if ~isstruct(block) || ~isfield(block, 'predictors'), continue; end

    preds = block.predictors;
    if isempty(preds), continue; end

    for j = 1:numel(preds)
        pname = preds{j}.name;
        varNames{end+1} = sprintf('%s_%s', comp, pname); %#ok<AGROW>
    end
end

% Make a 1xN table with those VariableNames (values not needed for labeling)
param_tbl = array2table(nan(1, numel(varNames)), 'VariableNames', varNames);

end


% ---------- your existing extractor (works unchanged) ----------
function [suffixes, prefixes] = extract_dep(results)

params = results.Properties.VariableNames;
parts_str = cellfun(@(s) split(s, '_'), params, 'UniformOutput', false);
suffixes = cellfun(@(parts) parts{end}, parts_str, 'UniformOutput', false);

% join everything except last back into one string (safer than a cell array)
prefixes = cellfun(@(parts) strjoin(parts(1:end-1), '_'), parts_str, 'UniformOutput', false);

end
