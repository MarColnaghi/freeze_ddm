function [xtrue, lbl, mask] = get_ground_truth_vector(model)

    components = fieldnames(model);
    xtrue_parts = cell(1, numel(components));
    base_names = {'sm', 'smp', 'fs', 'ln', 'ls', 'intercept'};

    for i = 1:numel(components)
        comp = model.(components{i});
        vec = nan(1, numel(base_names));
        zero_vec = zeros(1, numel(base_names));

        for idx_predictors = 1:size(comp.predictors, 2)
            zero_vec = strcmp(comp.predictors{idx_predictors}.name, base_names);
            vec(logical(zero_vec)) = comp.ground_truth(idx_predictors);
        end

        comp.mask = ~isnan(vec);
        xtrue_parts{i} = vec;
        xtrue_labels{i} = strcat(components{i}, "_", base_names);
        xtrue_mask{i} = comp.mask;
    end

    xtrue = horzcat(xtrue_parts{:});
    lbl = horzcat(xtrue_labels{:});
    mask = horzcat(xtrue_mask{:});
end
