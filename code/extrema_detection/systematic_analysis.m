threshold_imm = 3; threshold_mob = 3; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'extrema_detection/systematic_analysis', 'bouts_id', id_code);

sim_params.seed = 1;
sim_params.dt = 1/60;
sim_params.T = 30;
sim_params.time_vector = sim_params.dt:sim_params.dt:sim_params.T;
sim_params.z = 0;
sim_params.snr = 60;
sim_params.sigma_ed = 1;
sim_params.sigma_ac = 1;

% DDM params
sim_params.ddm.mu1 = 5;
sim_params.ddm.theta1 = 6;
sim_params.ddm.ndt = 0;

% Simulation settings
sim_params.kde_grid = 0:1/600:120;
sim_params.num_sims = 20000;

% Censoring/Truncation
sim_params.points.truncation = [];
sim_params.points.censoring = sim_params.T;

create_output_dirs(paths)

for idx_mu = [1, 3, 5, 10]
    for idx_theta = [2, 4, 6]
        for idx_snr = [30, 60]
       
            sim_params.ddm.mu1 = idx_mu;
            sim_params.ddm.theta1 = idx_theta;
            sim_params.ddm.ndt = 0;
            sim_params.snr = idx_snr;
            simulate_models('sim_params', sim_params, 'plot', true, 'store', true, 'paths', paths)

        end
    end
end

function create_output_dirs(paths)
% Ensure base directories exist
if ~exist(paths.fig, 'dir'), mkdir(paths.fig); end
if ~exist(paths.results, 'dir'), mkdir(paths.results); end

% Auto-incrementing run folder inside results
run_folders = dir(fullfile(paths.results, 'run*'));
run_nums = [];

for i = 1:length(run_folders)
    if run_folders(i).isdir
        tokens = regexp(run_folders(i).name, '^run(\d+)$', 'tokens');
        if ~isempty(tokens)
            run_nums(end+1) = str2double(tokens{1}{1}); %#ok<AGROW>
        end
    end
end

if isempty(run_nums)
    next_run = 1;
else
    next_run = max(run_nums) + 1;
end

run_name = sprintf('run%02d', next_run);
paths.results = fullfile(paths.results, run_name);
mkdir(paths.results);

% Also update figure path to match the new run
paths.fig = fullfile(paths.fig, run_name);
mkdir(paths.fig);

% Assign the updated paths back to base workspace (if needed)
assignin('caller', 'paths', paths);
end
