function pr = param_res(true_bound, points)

% -- Likelihood / Signal Parameters (Coarse) --
pr.dt = 1/60;       % 1ms for signal and PDE
pr.dx = 0.01;
pr.sigma_sq = 1.0;
%     pr.bound = 1.0;
pr.x0 = 0.0;

% PDE Grid Setup
x_min = -7;
grid_size = round((true_bound - x_min) / pr.dx) + 1;
start_idx = round((pr.x0 - x_min) / pr.dx) + 1;

% Save for objective function
pr.x_min = x_min;
pr.grid_size = grid_size;
pr.start_idx = start_idx;

% Truncation and T_max
pr.T_trunc = points.truncation;
pr.T_max = points.censoring - pr.dt;