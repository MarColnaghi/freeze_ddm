function [logL, kde] = kde_logL(rts, obs_rt, kde_grid)
    kde = ksdensity(rts, kde_grid, 'Function', 'pdf', 'Support','positive', 'Bandwidth', 1/60);
    p = interp1(kde_grid, kde, obs_rt, 'pchip');
    p = max(real(p), 1e-10);  % avoid log(0)
    logL = log(p);
end
