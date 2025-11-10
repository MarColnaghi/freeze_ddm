function out = evaluate_model(model, x, y)
% ------------------------------------------------------------
% y      : table with all predictors as columns
%         - scalar predictors: N×1 double
%         - time-series predictors: N×1 cell, each cell 1×T double
% x      : table with parameters (columns are coefficients)
% model  : struct with fields like model.mu1, model.theta1 ...
%          each containing:
%              .predictors {cell array of structs with .name}
%              .link       handle, e.g. @(z) z  or  @(z) exp(z)
% out    : table with same fields as `model`, each an [N x 1] vector
% ------------------------------------------------------------

    y_vars = y.Properties.VariableNames;
    x_vars = x.Properties.VariableNames;

    N = height(y);
    out = table();

    param_names = fieldnames(model);

    for i = 1:numel(param_names)
        pname = param_names{i};
        block = model.(pname);

        % All coefficients belonging to this parameter block
        tf_block = startsWith(x_vars, pname);
        x_block  = x(:, tf_block);   % table of coefficients for this block

        % Initialize linear predictor
        eta = zeros(N, 1);

        % Loop over predictors in this block
        for j = 1:numel(block.predictors)
            pred_name = block.predictors{j}.name;

            % Column of y corresponding to this predictor
            if ~ismember(pred_name, y_vars)
                error('Predictor "%s" not found in y.', pred_name);
            end
            Xcol = y.(pred_name);   % could be numeric or cell

            % Coefficient(s) corresponding to this predictor in x_block
            tf_beta = endsWith(x_block.Properties.VariableNames, pred_name);
            if ~any(tf_beta)
                error('Coefficient for predictor "%s" not found in x.', pred_name);
            end
            beta = table2array(x_block(:, tf_beta));  % row or scalar

            % ---- Handle scalar vs time-series predictors ----
            if iscell(Xcol)
                % Time-series predictor: N×1 cell, each 1×T double
                Xi = cell2mat(Xcol);  % N×T

                if isscalar(beta)
                    % same weight for all time points
                    eta = eta + Xi * beta;          % N×1
                else
                    % one weight per time point: beta should be 1×T or T×1
                    beta = beta(:);                % T×1
                    eta = eta + Xi * beta;         % N×1: row-wise dot product
                end

            else
                % Ordinary numeric predictor: assume N×1
                if isvector(beta)
                    beta = beta(1);   % safety: take scalar if they gave a vector
                end
                eta = eta + Xcol * beta;  % N×1
            end
        end

        % Apply link and store in output table
        out.(pname) = block.link(eta);
    end
end
