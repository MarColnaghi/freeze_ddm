%% The Trapezoidal Cure: Scaling for High Dimensions
clear; clc; close all;

n_sim = 10000;
k = 5; % High number of predictors

% --- Helper: Function to sample from a Trapezoid centered at 0 ---
% A trapezoid is the sum of two uniform distributions: U(-a, a) + U(-b, b)
draw_trap = @(n, k, scale) (-scale + (2*scale)*rand(n,k)) + ...
                            (-scale/2 + scale*rand(n,k));

% --- VERSION A: Unscaled Trapezoid (The "Sick" Model) ---
% The width is fixed at 5, regardless of how many variables we have.
betas_wide = draw_trap(n_sim, k, 5); 
z_wide = sum(betas_wide, 2);
prob_wide = 1 ./ (1 + exp(-z_wide));

% --- VERSION B: Scaled Trapezoid (The "Cure") ---
% We shrink the width by 1/sqrt(k) to keep the output sane.
scale_factor = 5 / sqrt(k); 
betas_scaled = draw_trap(n_sim, k, scale_factor);
z_scaled = sum(betas_scaled, 2);
prob_scaled = 1 ./ (1 + exp(-z_scaled));

% --- Visualization ---
figure('Color', 'w', 'Position', [100, 100, 1000, 400]);

subplot(1,2,1);
histogram(prob_wide, 30, 'FaceColor', [0.8 0.3 0.3], 'EdgeColor', 'w');
title(['Wide Trapezoid (k=', num2str(k), ')']);
xlabel('Predicted Probability'); ylabel('Simulations');
subtitle('Problem: Polarized at 0 and 1');

subplot(1,2,2);
histogram(prob_scaled, 30, 'FaceColor', [0.3 0.8 0.3], 'EdgeColor', 'w');
title(['Scaled Trapezoid (k=', num2str(k), ')']);
xlabel('Predicted Probability');
subtitle('Cure: Balanced "Maybe" Predictions');

sgtitle('Prior Predictive Check: Trapezoidal Priors');


%%
%% The "Deep Cure" for Trapezoidal Bimodality
clear; clc; close all;

n_sim = 10000;
k = 5; % Number of variables

% Function to draw from a trapezoid (Sum of two Uniforms)
% Width 'a' controls the flat top, 'b' controls the total spread
draw_trap = @(n, k, a, b) (-a + 2*a*rand(n,k)) + (-b + 2*b*rand(n,k));

% --- EXPERIMENT 1: The "Weak" Scale (Causes Bimodality) ---
% Even if we scale by 1/sqrt(k), if the starting width is too high, 
% the sum is too wide for the Sigmoid function.
scale_weak = 2 / sqrt(k); 
z_weak = sum(draw_trap(n_sim, k, scale_weak, scale_weak/2), 2);
prob_weak = 1 ./ (1 + exp(-z_weak));

% --- EXPERIMENT 2: The "Outcome-Driven" Scale (The Real Cure) ---
% We want the total sum (z) to have a Standard Deviation of roughly 1 or 2.
% We scale the individual priors much more aggressively.
scale_strong = 0.5 / sqrt(k); 
z_strong = sum(draw_trap(n_sim, k, scale_strong, scale_strong/2), 2);
prob_strong = 1 ./ (1 + exp(-z_strong));

% --- Visualization ---
figure('Color', 'w', 'Position', [100, 100, 1000, 450]);

subplot(1,2,1);
histogram(prob_weak, 30, 'FaceColor', [0.8 0.4 0.4]);
title('Still Bimodal (Scale too wide)');
xlabel('Predicted Probability'); ylabel('Count');
subtitle('The sum of priors is still "overpowering" the model');

subplot(1,2,2);
histogram(prob_strong, 30, 'FaceColor', [0.4 0.7 0.4]);
title('Sensible/Unimodal (Strongly Scaled)');
xlabel('Predicted Probability');
subtitle('The "Cure": Forcing the outcome to stay near 0.5');

sgtitle('Prior Predictive Check: Adjusting the "Volume" of the Priors');