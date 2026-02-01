function [lb, plb, pub, ub, prior_info] = extract_bounds_from_model(model)

lb = [];
ub = [];
prior_info = struct('type', {}, 'lambda', {});

param_blocks = fieldnames(model);
base_names = {'sm', 'smp', 'fs', 'ln', 'ls', 'tsl', 'intercept'};

for i = 1:numel(param_blocks)
    block = model.(param_blocks{i});
    preds = block.predictors;

    block_lb = nan(1, numel(base_names));
    block_ub = nan(1, numel(base_names));
    block_prior = repmat(struct('type','trapezoidal','lambda',[]), ...
                         1, numel(base_names));

    for j = 1:numel(preds)
        pred = preds{j};
        pred_name = pred.name;
        bounds = pred.bounds;

        idx = find(strcmp(base_names, pred_name));
        if isempty(idx)
            continue
        end

        block_lb(idx) = bounds(1);
        block_ub(idx) = bounds(2);

        % ---- NEW: read prior info if present
        if isfield(pred, 'prior')
            block_prior(idx).type = pred.prior;

            if strcmp(pred.prior, 'exponential')
                block_prior(idx).lambda = pred.lambda;
                
            end
        end
    end

    % Append only active parameters
    active = ~isnan(block_lb);

    lb = [lb, block_lb(active)];
    ub = [ub, block_ub(active)];
    prior_info = [prior_info, block_prior(active)];
end

% Plausible bounds
plb = lb + 0.1;
pub = ub - 0.1;

fixed_idx = (lb == ub);
plb(fixed_idx) = lb(fixed_idx);
pub(fixed_idx) = ub(fixed_idx);


end
