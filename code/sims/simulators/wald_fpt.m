function [pdf, cdf] = wald_fpt(t, a, mu, sigma)
% WALD_FPT  Exact first-passage density/CDF of Brownian motion with drift to a
%           SINGLE absorbing barrier -- the Wald (inverse Gaussian) distribution.
%
%   [pdf, cdf] = wald_fpt(t, a, mu, sigma)
%
%   t     : times (vector, > 0)
%   a     : barrier distance from the start point (a = bound - x0), a > 0
%   mu    : constant drift RATE
%   sigma : diffusion SD (default 1)
%
%   For  dx = mu dt + sigma dW,  x(0) = 0,  absorbed at x = a > 0, NO lower
%   boundary, the first-passage time T has density
%
%       f(t) = a / (sigma * sqrt(2*pi*t^3)) * exp( -(a - mu*t)^2 / (2*sigma^2*t) )
%
%   and distribution
%
%       F(t) = Phi( (mu*t - a)/(sigma*sqrt(t)) )
%              + exp(2*mu*a/sigma^2) * Phi( -(mu*t + a)/(sigma*sqrt(t)) )
%
%   This is EXACT -- no grid, no sampling -- so it is ground truth for the
%   Crank-Nicolson solver (see fpe_cn.m) whenever the drift is constant and there
%   is no leak. Two conditions must hold for the comparison to be fair:
%     * lambda = 0 (the Wald has no leak; with leak the process is
%       Ornstein-Uhlenbeck and has no such closed form), and
%     * the solver's reflecting floor x_min must sit far enough below that
%       essentially no mass reaches it -- the Wald assumes NO lower boundary.
%
%   For mu <= 0 the barrier is not hit with probability 1, and F(inf) =
%   exp(2*mu*a/sigma^2) < 1; the formulas remain valid (defective distribution).
%
%   Matches bayes_fpe's wald_log_pdf / wald_cdf (tv.py:181-193), which are the
%   sigma = 1 special case.

    if nargin < 4 || isempty(sigma), sigma = 1; end
    t = t(:);
    pdf = zeros(size(t));
    cdf = zeros(size(t));

    ok = t > 0;
    tt = t(ok);

    pdf(ok) = a ./ (sigma * sqrt(2*pi*tt.^3)) .* ...
              exp( -(a - mu*tt).^2 ./ (2 * sigma^2 * tt) );

    if nargout > 1
        z1 = (mu*tt - a) ./ (sigma * sqrt(tt));
        z2 = -(mu*tt + a) ./ (sigma * sqrt(tt));
        % exp(2*mu*a/sigma^2) can overflow for large mu*a; pair it with the log
        % of the normal CDF so the product stays finite.
        cdf(ok) = normcdf(z1) + exp( 2*mu*a/sigma^2 + log(max(normcdf(z2), realmin)) );
        cdf = min(max(cdf, 0), 1);
    end
end
