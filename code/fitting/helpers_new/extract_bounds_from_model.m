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

    for i = 1:numel(param_blocks)
        block = model.(param_blocks{i});
        preds = block.predictors;

        % Only keep bounds of active predictors

        for j = 1:numel(preds)
            bounds = preds{j}.bounds;
            lb(1, end + 1) = bounds(1);
            ub(1, end + 1) = bounds(2);
        end
    end
plb = lb + 0.1;
pub = ub - 0.1;
end

