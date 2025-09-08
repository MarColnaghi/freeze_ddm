function [model] = generate_masks(model)

param_blocks = fieldnames(model);  % e.g. {'mu1', 'theta1', 'tndt1'}

for i = 1:numel(param_blocks)
    block = model.(param_blocks{i});
    preds = block.predictors;
    mask = zeros(1,6);

    for idx_preds = 1:length(preds)
        if strcmp(preds{idx_preds}.name, 'sm')
            mask(1) = 1;
        elseif strcmp(preds{idx_preds}.name, 'fs')
            mask(2) = 1;
        elseif strcmp(preds{idx_preds}.name, 'ln')
            mask(3) = 1;
        elseif strcmp(preds{idx_preds}.name, 'ls')
            mask(4) = 1;
        elseif strcmp(preds{idx_preds}.name, 'intercept')
            mask(5) = 1;
        elseif strcmp(preds{idx_preds}.name, 'smp')
            mask(6) = 1;
        end

    end
    model.(param_blocks{i}).mask = mask;

end