function [logL, f_kde_scaled] = kde_logL_censored(rts, obs_rt, kde_grid)
    
    is_censored = isnan(rts);
    is_observed = ~is_censored;
    rts_obs = rts(is_observed);      % Uncensored RTs
    n_total = numel(rts);
    n_obs = sum(is_observed);
    n_cens = sum(is_censored);
    
    kde = ksdensity(rts_obs, kde_grid, 'Function', 'pdf', 'Bandwidth', 1/30);
    dx = kde_grid(2) - kde_grid(1);
    area_kde = sum(kde * dx);  % Area under the KDE curve

    p_obs = n_obs / n_total;      % Probability mass below tcensor

    f_kde_scaled = kde * p_obs / area_kde;

    % Plot scaled KDE for observed RTs
    f_kde_scaled(kde_grid > max(rts)) = 1 - trapz(kde_grid, f_kde_scaled);
    
%         figure;
%     hold on;
%     plot(kde_grid, kde, 'b', 'LineWidth', 2);
%     plot(kde_grid, f_kde_scaled, 'r', 'LineWidth', 2);
% 
%     histogram(rts, 1/120:1/10:12,'FaceColor', 'r', 'Normalization', 'pdf')
%     histogram(rts(~is_censored), 1/120:1/10:12, 'FaceColor', 'b', 'Normalization', 'pdf')
% 
%     % Plot censored point mass at tcensor
%     stem(10.5, p_cens / dx, 'r', 'LineWidth', 2);  % scaled as a spike
% 
%     xlabel('Reaction Time');
%     ylabel('Estimated PDF');
%     title('Mixed KDE of DDM with Right-Censoring');

    p = interp1(kde_grid, f_kde_scaled, obs_rt, 'pchip');
    p = max(real(p), 1e-6);  % avoid log(0)
    logL = log(p);
end
