function plot_model_performance(dataA, idsA, dataB, idsB, M)
    % Use the true parameters for a "clean" comparison
    p = 0.7; mu = [0, 5]; sig = [1, 1];
    
    llA_on_A = zeros(M, 1);
    llB_on_A = zeros(M, 1);
    llA_on_B = zeros(M, 1);
    llB_on_B = zeros(M, 1);

    for j = 1:M
        % Extract data for subject j
        xA = dataA(idsA == j);
        xB = dataB(idsB == j);
        
        % Model A Subject Scores (Independent)
        llA_on_A(j) = sum(log(p*normpdf(xA, mu(1), sig(1)) + (1-p)*normpdf(xA, mu(2), sig(2))));
        llA_on_B(j) = sum(log(p*normpdf(xB, mu(1), sig(1)) + (1-p)*normpdf(xB, mu(2), sig(2))));
        
        % Model B Subject Scores (Loyalist)
        % Note: We use the stable Log-Sum-Exp logic here
        lB1_A = log(p) + sum(log(normpdf(xA, mu(1), sig(1))));
        lB2_A = log(1-p) + sum(log(normpdf(xA, mu(2), sig(2))));
        llB_on_A(j) = max(lB1_A, lB2_A) + log(exp(lB1_A-max(lB1_A,lB2_A)) + exp(lB2_A-max(lB1_A,lB2_A)));
        
        lB1_B = log(p) + sum(log(normpdf(xB, mu(1), sig(1))));
        lB2_B = log(1-p) + sum(log(normpdf(xB, mu(2), sig(2))));
        llB_on_B(j) = max(lB1_B, lB2_B) + log(exp(lB1_B-max(lB1_B,lB2_B)) + exp(lB2_B-max(lB1_B,lB2_B)));
    end

    %% Plotting
    figure('Color', 'w', 'Position', [100, 100, 1100, 500]);
    
    subplot(1,2,1);
    hold on;
    plot(1:M, llA_on_A, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Model A (Shifter)');
    plot(1:M, llB_on_A, 'bx', 'LineWidth', 2, 'DisplayName', 'Model B (Loyalist)');
    title('Fitting to SHIFTER Data');
    ylabel('Log-Likelihood per Subject'); xlabel('Subject ID');
    legend('Location', 'southoutside', 'Orientation', 'horizontal'); grid on;

    subplot(1,2,2);
    hold on;
    plot(1:M, llA_on_B, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Model A (Shifter)');
    plot(1:M, llB_on_B, 'bx', 'LineWidth', 2, 'DisplayName', 'Model B (Loyalist)');
    title('Fitting to LOYALIST Data');
    ylabel('Log-Likelihood per Subject'); xlabel('Subject ID');
    legend('Location', 'southoutside', 'Orientation', 'horizontal'); grid on;
end