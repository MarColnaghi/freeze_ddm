function [lb, plb, pub, ub, prior_info] = extract_bounds_from_model(model)

lb = [];
ub = [];
plb = [];
pub = [];
prior_info = struct('type', {}, 'lambda', {});

param_blocks = fieldnames(model);
base_names = {'sm', 'smp', 'fs', 'ln', 'ls', 'tsl', 'intercept'};

for i = 1:numel(param_blocks)
    block = model.(param_blocks{i});
    preds = block.predictors;

    block_lb = nan(1, numel(base_names));
    block_ub = nan(1, numel(base_names));

    block_plb = nan(1, numel(base_names));
    block_pub = nan(1, numel(base_names));

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

        if numel(bounds) == 2

            block_lb(idx) = bounds(1);
            block_ub(idx) = bounds(2);

        elseif numel(bounds) == 4

            block_lb(idx) = bounds(1);
            block_ub(idx) = bounds(4);

            block_plb(idx)  = bounds(2);
            block_pub(idx) = bounds(3);

        end

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

    plb = [plb, block_plb(active)];
    pub = [pub, block_pub(active)];

    prior_info = [prior_info, block_prior(active)];
end

% We rewrite the pub and plb in the case those are not specified.

if numel(bounds) == 2

    % Plausible bounds
    plb = lb + 0.1;
    pub = ub - 0.1;

end

fixed_idx = (lb == ub);
plb(fixed_idx) = lb(fixed_idx);
pub(fixed_idx) = ub(fixed_idx);


end
