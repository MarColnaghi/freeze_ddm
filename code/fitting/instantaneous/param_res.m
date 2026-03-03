function pr = param_res(true_bound, varargin)

opt = inputParser;
addParameter(opt, 'points', struct('censoring', 10.5, 'truncation', 0));
addParameter(opt, 'dt', 1/60);
parse(opt, varargin{:});

% -- Likelihood / Signal Parameters (Coarse) --
pr.dt = opt.Results.dt;       % 1ms for signal and PDE
pr.dx = 0.001;
pr.sigma_sq = 1.0;
%     pr.bound = 1.0;
pr.x0 = 0.0;

% PDE Grid Setup
x_min = -8;
grid_size = round((true_bound - x_min) / pr.dx) + 1;
start_idx = round((pr.x0 - x_min) / pr.dx) + 1;

% Save for objective function
pr.x_min = x_min;
pr.grid_size = grid_size;
pr.start_idx = start_idx;

% Truncation and T_max
pr.T_trunc = opt.Results.points.truncation;
pr.T_max = opt.Results.points.censoring;