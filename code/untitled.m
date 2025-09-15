%% outlier_model_demo.m
% Outlier detection model: builds decision-time PMF from hazard,
% shifted log-logistic delay PDF/CDF, and response-time PDF/CDF via convolution.

clear; clc;

%% ------------------------ Parameters & stimulus -------------------------
fs = 60;                         % Hz (video frames)
T  = 2;                          % total duration (s)
t  = (0:1/fs:T)';                % time grid (column)

% Example stimulus s_i (log TF); replace with your actual stimulus trace
s = 1.0*sin(2*pi*2*t);     % small modulation

% Sensory / decision params
sigma = 0.25;                    % Gaussian noise std
b     = 0.8;                     % decision bound

% Shifted log-logistic non-decision delay params (alpha>0, beta>0, gamma>1)
alpha = 0.08; beta = 0.18; gamma = 2.2;

% Response-time grid (can be denser than t if you want smoother curves)
r_grid = t;                      % here we just reuse t

%% ---------------------- Compute distributions --------------------------
res = outlier_model_response_distribution(s, t, sigma, b, alpha, beta, gamma, r_grid);

fprintf('Mass of pD over finite decisions: %.6f\n', sum(res.pD));
fprintf('Final survival S(end): %.6f (mass at D=∞ = S(end))\n', res.survival(end));
fprintf('Final FR(end): %.6f (some mass may remain for D=∞)\n', res.FR(end));

%% ----------------------------- Plots -----------------------------------
figure; plot(res.t, res.pD, '-','LineWidth',1.5);
xlabel('time (s)'); ylabel('p_D(t_i)'); title('Decision-time PMF'); grid on;

figure; plot(res.r, res.pR, '-','LineWidth',1.5);
xlabel('time (s)'); ylabel('p_R(r)');  title('Response-time PDF'); grid on;

figure; plot(res.r, res.FR, '-','LineWidth',1.5);
xlabel('time (s)'); ylabel('F_R(r)');  title('Response-time CDF'); grid on;

%% ========================= Local functions =============================
function out = outlier_model_response_distribution(s, t, sigma, b, alpha, beta, gamma, r_grid)
% Wrapper: builds pD from hazard and then pR, FR on r_grid.
    [pD, h, S] = decision_pmf_from_hazard(s, t, sigma, b);
    if nargin < 8 || isempty(r_grid), r_grid = t; end
    [pR, FR] = response_pdf_cdf(t, pD, r_grid, alpha, beta, gamma);
    out = struct('t', t, 'pD', pD, 'hazard', h, 'survival', S, ...
                 'r', r_grid, 'pR', pR, 'FR', FR);
end

function [pD, h, S] = decision_pmf_from_hazard(s, t, sigma, b)
% Decision-time PMF from discrete hazard.
% h_i = 1 - Phi((b - s_i)/sigma)
% pD_i = h_i * prod_{j<i} (1 - h_j)
% S_i  = prod_{j<=i} (1 - h_j)
    s = s(:); t = t(:); %#ok<NASGU> (t included for completeness/signature)
    z = (b - s)./sigma;
    Phi = 0.5*(1 + erf(z./sqrt(2)));   % standard normal CDF
    h   = 1 - Phi;                     % hazard per frame
    % Numerical stability via log-survival accumulation
    log1mh   = log(max(1 - h, realmin));
    logSprev = [0; cumsum(log1mh(1:end-1))];  % log S_{i-1}
    pD       = h .* exp(logSprev);
    S        = exp(cumsum(log1mh));           % S_i
end

function [pR, FR] = response_pdf_cdf(t, pD, r_grid, alpha, beta, gamma)
% p_R(r_k) = sum_{i: t_i <= r_k} pD(t_i) * fdelta(r_k - t_i)
% F_R(r_k) = sum_{i: t_i <= r_k} pD(t_i) * Fdelta(r_k - t_i)
    t      = t(:);
    pD     = pD(:);
    r_grid = r_grid(:);

    nT = numel(t);
    nR = numel(r_grid);

    pR = zeros(nR,1);
    FR = zeros(nR,1);

    for i = 1:nT
        di = t(i);
        if pD(i) == 0, continue; end
        mask = r_grid >= di;                % only r >= d contribute
        if ~any(mask), continue; end
        tau = r_grid(mask) - di;            % delays r - d
        [fdelta, Fdelta] = delay_loglogistic_pdf_cdf(tau, alpha, beta, gamma);
        pR(mask) = pR(mask) + pD(i) * fdelta;
        FR(mask) = FR(mask) + pD(i) * Fdelta;
    end
end

function [f, F] = delay_loglogistic_pdf_cdf(delta, alpha, beta, gamma)
% Shifted log-logistic PDF/CDF (zero for delta < alpha).
% x = (delta - alpha)/beta
% f(delta) = (gamma/beta) * x^(gamma-1) / (1 + x^gamma)^2
% F(delta) = 1 / (1 + x^{-gamma})
    delta = delta(:);
    f = zeros(size(delta));
    F = zeros(size(delta));
    x = (delta - alpha) ./ beta;
    valid = x >= 0;
    xv = x(valid);
    denom = (1 + xv.^gamma);
    f(valid) = (gamma./beta) .* (xv.^(gamma-1)) ./ (denom.^2);
    F(valid) = 1 ./ denom;
end

function delta = delay_loglogistic_icdf(u, alpha, beta, gamma)
% Inverse CDF for sampling delays (optional utility).
% Delta = alpha + beta * (u/(1-u))^(1/gamma)
    delta = alpha + beta .* (u./(1-u)).^(1./gamma);
end
