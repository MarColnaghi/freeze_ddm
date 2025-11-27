%% --- 4. Run Monte Carlo Simulation (C++ MEX) ---

if do_sims
    
N_trials = 100000;
fprintf('Running %d Monte Carlo simulation trials (C++)...\n', N_trials);


% Signal
 %% --- 2. Generate the Drift Time Series ---
fprintf('Generating drift time series...\n');
% (Using the same stochastic generation code as before)
std_drift = 0.2; 
correlation_time_ms = 100; 
peak_amplitude = 2.0; 
peak_probability = 0.001; 
peak_decay_ms = 10; 
sampling_interval_ms = 10;
filter_length = round(correlation_time_ms / 1000 / params.dt); 
if filter_length < 1, filter_length = 1; end
filter_coeffs = ones(1, filter_length) / filter_length; 
correlated_noise = filter(filter_coeffs, 1, randn(1, time_steps)); 
correlated_noise = (correlated_noise - mean(correlated_noise)) / std(correlated_noise); 
base_process = params.mean_drift + std_drift * correlated_noise; 
decay_constant = round(peak_decay_ms / 1000 / params.dt); 
if decay_constant < 1, decay_constant = 1; end 
peak_triggers = find(rand(1, time_steps) < peak_probability); 
peak_train = zeros(1, time_steps); 
for trigger_idx = peak_triggers 
    decay_length = min(5*decay_constant, time_steps - trigger_idx + 1); 
    decay_shape = exp(-(0:decay_length-1) / decay_constant); 
    peak_train(trigger_idx : trigger_idx + decay_length - 1) = peak_train(trigger_idx : trigger_idx + decay_length - 1) + decay_shape; 
end
if ~isempty(peak_triggers) && max(peak_train) > 0, peak_train = peak_train / max(peak_train) * peak_amplitude; end
drifts_continuous = base_process + peak_train; 
sampling_interval_steps = round(sampling_interval_ms / 1000 / params.dt); 
if sampling_interval_steps < 1, sampling_interval_steps = 1; end
drifts = zeros(1, time_steps); 
for start_idx = 1 : sampling_interval_steps : time_steps
    measured_drift = drifts_continuous(start_idx); 
    end_idx = min(start_idx + sampling_interval_steps - 1, time_steps); 
    drifts(start_idx : end_idx) = measured_drift; 
end
fprintf('Drift generation complete.\n\n');

% drifts(:) =  params.mean_drift;

% Prepare parameters [dt, sigma_w_sq, x0, theta, lambda]
dt_for_sim = .00005;
time_vec_for_sim = 0:dt_for_sim:params.T_max;
drifts_for_sim = interp1(time_vec,drifts,time_vec_for_sim,'nearest');

sim_params = [dt_for_sim, params.sigma_w_sq, params.x0, params.theta, params.lambda];

tic;

% --- CALL THE MEX FUNCTION ---
fpt_results = sim_ddm_parallel(drifts_for_sim, sim_params, N_trials);
% -----------------------------

% Post-processing (Handle non-crossings)
% If it returns NaN, it means it didn't cross. Set to T_max.
fpt_results(isnan(fpt_results)) = time_vec(end);

elapsed_time = toc;
fprintf('Simulation finished in %.4f seconds.\n\n', elapsed_time);

end