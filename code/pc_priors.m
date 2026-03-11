%% Step 1: PC Prior Setup (Principles 3 & 4)
% Define the scaling statement: Prob(p > U) = alpha [cite: 328]
U = 0.2;      % Upper bound for mixture probability p [cite: 329]
alpha = 0.1;  % Weight put on this tail event [cite: 329]

% Calculate Lambda (Decay Rate) [cite: 359]
% PC Prior is a truncated exponential on the distance scale d [cite: 318, 1083]
% Here we assume d is roughly linear with p for the mixture weight [cite: 916]
lambda_pc = -log(alpha) / U; 

% PC Prior Density Function [cite: 1083]
pc_prior_pdf = @(p) (lambda_pc * exp(-lambda_pc * p)) / (1 - exp(-lambda_pc));

%% Step 2: Generate Synthetic Data
% We use the Inverse Gaussian distribution as a proxy for DDM First Passage Time
num_samples = 400;

% Process 1 (Fast/Short Bouts): drift=2.0, boundary=1.0
mu1 = 1.0 / 2.0; lambda1 = 1.0^2; 
% Process 2 (Slow/Long Bouts): drift=0.5, boundary=1.0
mu2 = 1.0 / 0.5; lambda2 = 1.0^2;

% Scenario A: Data is 100% Process 1 (The Base Model [cite: 190])
data_base = random('InverseGaussian', mu1, lambda1, [num_samples, 1]);

% Scenario B: Data is a 30/70 Mixture (The Flexible Extension [cite: 193])
idx = rand(num_samples, 1) > 0.6;
data_mix = zeros(num_samples, 1);
data_mix(idx) = random('InverseGaussian', mu1, lambda1, [sum(idx), 1]);
data_mix(~idx) = random('InverseGaussian', mu2, lambda2, [sum(~idx), 1]);

%% Step 3: Compute Likelihoods and Posteriors
p_grid = linspace(0.001, 0.999, 100);
[post_base, post_mix] = deal(zeros(size(p_grid)));

for i = 1:length(p_grid)
    p = p_grid(i);
    
    % Likelihood for Scenario A
    L_base = (1-p)*pdf('InverseGaussian', data_base, mu1, lambda1) + ...
               p*pdf('InverseGaussian', data_base, mu2, lambda2);
    post_base(i) = sum(log(L_base + eps)) + log(pc_prior_pdf(p));
    
    % Likelihood for Scenario B
    L_mix = (1-p)*pdf('InverseGaussian', data_mix, mu1, lambda1) + ...
              p*pdf('InverseGaussian', data_mix, mu2, lambda2);
    post_mix(i) = sum(log(L_mix + eps)) + log(pc_prior_pdf(p));
end

% Normalize Posteriors
post_base = exp(post_base - max(post_base)); post_base = post_base / sum(post_base);
post_mix = exp(post_mix - max(post_mix)); post_mix = post_mix / sum(post_mix);

%% Step 4: Visualization
figure('Color', 'w', 'Position', [100 100 800 600]);

% Plot Prior
subplot(2,1,1);
plot(p_grid, pc_prior_pdf(p_grid), 'r--', 'LineWidth', 2);
title('PC Prior: Penalising Model Complexity (Base Model p=0) [cite: 297]');
ylabel('Density'); xlabel('Mixture Probability (p)'); grid on;

% Plot Results
subplot(2,1,2); hold on;
plot(p_grid, post_base, 'b', 'LineWidth', 2, 'DisplayName', 'Posterior (True p=0)');
plot(p_grid, post_mix, 'k', 'LineWidth', 2, 'DisplayName', 'Posterior (True p=0.3)');
legend('Location', 'northeast');
title('Posterior Comparison: Base Model vs. Mixture');
xlabel('Mixture Probability (p)'); ylabel('Normalized Density'); grid on;