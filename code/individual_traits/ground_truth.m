function outlier_fragility_demo()
    
% Parameters for our two Gaussians
    mu1 = 0; sig1 = 1;
    mu2 = 5; sig2 = 1;
    p = 0.5; % Equal prior for both groups
    
    % Generate a "mostly clean" subject
    % 19 points are near Mean 1 (0)
    n_clean = 10;
    clean_data = mu1 + sig1 * randn(n_clean, 1);
    
    % Range of outlier valu
    % es for the 20th point
    outlier_range = linspace(-5, 15, 100);
    
    ll_A = zeros(size(outlier_range));
    ll_B = zeros(size(outlier_range));
    
    for k = 1:length(outlier_range)
        current_data = [clean_data; outlier_range(k)];
        
        % --- Model A (Shifter) ---
        % Sum of log-mixtures: Each point is handled independently
        % log(p*N1 + (1-p)*N2) for every point
        probs_A = p * normpdf(current_data, mu1, sig1) + ...
                  (1-p) * normpdf(current_data, mu2, sig2);
        ll_A(k) = sum(log(probs_A));
        
        % --- Model B (Loyalist) ---
        % Log of the sum of products: The whole block is one type
        % log(p * prod(N1) + (1-p) * prod(N2))
        prob_block1 = p * prod(normpdf(current_data, mu1, sig1));
        prob_block2 = (1-p) * prod(normpdf(current_data, mu2, sig2));
        
        % We use a small epsilon to avoid log(0) for visualization
        ll_B(k) = log(max(prob_block1 + prob_block2, 1e-300));
    end
    
    %% Plotting the Fragility
    figure('Color', 'w', 'Position', [100 100 800 500]);
    plot(outlier_range, ll_A, 'r', 'LineWidth', 2); hold on;
    plot(outlier_range, ll_B, 'b', 'LineWidth', 2);
    
    % Add a vertical line where the outlier hits the other Gaussian
    line([5 5], get(gca,'YLim'), 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
    text(5.2, min(ll_A), 'Mean of Gaussian 2', 'Color', [0.5 0.5 0.5]);
    
    xlabel('Value of the 20th Sample (Outlier)');
    ylabel('Log-Likelihood for the Subject');
    title('Fragility: Model A (Shifter) vs Model B (Loyalist)');
    legend('Model A: Hedges bets (Independent)', 'Model B: Bets the house (Trait)');
    grid on;
end