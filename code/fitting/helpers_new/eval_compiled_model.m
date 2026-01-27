function out = eval_compiled_model(blocks, x)
% FAST NUMERIC EVALUATION
%
% x : parameter vector

    % Infer N once
    N = numel(blocks(1).X{1});
    out = struct;

    for i = 1:numel(blocks)
        eta = zeros(N,1);

        for j = 1:numel(blocks(i).X)
            beta = x(blocks(i).beta_inds(j));
            Xj   = blocks(i).X{j};

            if blocks(i).is_cell(j)
                % Time-series predictor: N x T
                Xi = cell2mat(Xj);
                eta = eta + Xi * beta;
            else
                % Scalar predictor
                eta = eta + Xj * beta;
            end
        end

        out.(blocks(i).name) = blocks(i).link(eta);
    end
end
