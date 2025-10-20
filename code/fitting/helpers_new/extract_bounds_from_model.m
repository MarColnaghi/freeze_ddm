function [lb, plb, pub, ub] = extract_bounds_from_model(model)
% ------------------------------------------------------------
% Extracts lower and upper bounds for all active predictors
% in the model, in the order they would appear in x
%
% Returns:
%   lb: [N x 1] vector of lower bounds
%   ub: [N x 1] vector of upper bounds
% ------------------------------------------------------------

lb = [];
ub = [];

param_blocks = fieldnames(model);  % e.g. {'mu1', 'theta1', 'tndt1'}
base_names = {'sm', 'smp', 'fs', 'ln', 'ls', 'intercept'};

for i = 1:numel(param_blocks)
    block = model.(param_blocks{i});
    preds = block.predictors;

    % Temporary storage for one parameter block
    block_lb = nan(1, numel(base_names));
    block_ub = nan(1, numel(base_names));

    % Match predictor names to base_names
    for j = 1:numel(preds)
        pred_name = preds{j}.name;   % assuming each predictor has a .name field like 'sm', 'fs', etc.
        bounds = preds{j}.bounds;

        % Find where this predictor fits in base_names
        idx = find(strcmp(base_names, pred_name));
        if ~isempty(idx)
            block_lb(idx) = bounds(1);
            block_ub(idx) = bounds(2);
        end
    end

    % Concatenate block bounds in the correct order
    lb = [lb, block_lb];
    ub = [ub, block_ub];
end
lb = lb(~isnan(lb));
ub = ub(~isnan(ub));
plb = lb + 0.1;
pub = ub - 0.1;


end
