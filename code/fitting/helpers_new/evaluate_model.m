function out = evaluate_model(model, x, y)
% ------------------------------------------------------------
% y      : table with all predictors as columns
% x      : vector of all parameters (in block order)
% model  : struct with fields like model.mu1, model.theta1 ...
%          each containing:
%              .predictors {cell array of structs with .name}
%              .link       handle, e.g. @(z) z  or  @(z) exp(z)
% out    : table with same fields as `model`, each an [N x 1] vector
% ------------------------------------------------------------

    y_vars   = y.Properties.VariableNames;
    x_vars   = x.Properties.VariableNames;
    X_full = y{:,:};                  % [N x P]
    out    = table();

    param_names = fieldnames(model);
    idx_start = 1;

    for i = 1:numel(param_names)
        pname = param_names{i};
        block = model.(pname);

        n_params = numel(block.predictors);           % number of coefficients in this block
        idx_end  = idx_start + n_params - 1;
        tf_start = startsWith(x_vars, pname);
        x_block  = x(:, tf_start);              % slice out relevant coefficients

        % Build full-length coefficient vector aligned to column order in y
        coeff = zeros(1, numel(y_vars));
        for j = 1:n_params
            pred_name = block.predictors{j}.name;
            idx = strcmp(y_vars, pred_name);
            tf_end = endsWith(x_block.Properties.VariableNames, pred_name);
            coeff(idx) = table2array(x_block(:, tf_end));
        end

        eta = X_full * coeff';
        out.(pname) = block.link(eta);

        idx_start = idx_end + 1;   % advance to next block
    end
end

