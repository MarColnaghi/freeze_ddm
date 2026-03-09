function main_demonstration()
    M = 30; n_j = 10;
    
    % 1. Generate Datasets
    [dataA, idsA, dataB, idsB] = generate_comparison_datasets(M, n_j);
    
    % 2. Run Cross-Fit Analysis
    results = cross_fit_demonstration(dataA, idsA, dataB, idsB, M);
    
    % 3. Extract and Plot Identities (The Identity Histogram)
    plot_identity_distributions(results, dataA, idsA, dataB, idsB, M);
end

function [dataA, idsA, dataB, idsB] = generate_comparison_datasets(M, n_j)
    p = 0.7; mu = [0, 5]; sig = [1, 1];
    
    % Process 1: Shifters
    dataA = zeros(M * n_j, 1); idsA = zeros(M * n_j, 1);
    for j = 1:M
        for i = 1:n_j
            type = (rand() > p) + 1;
            idx = (j-1)*n_j + i;
            dataA(idx) = mu(type) + sig(type) * randn();
            idsA(idx) = j;
        end
    end

    % Process 2: Loyalists
    dataB = zeros(M * n_j, 1); idsB = zeros(M * n_j, 1);
    for j = 1:M
        type = (rand() > p) + 1;
        idx = (j-1)*n_j + 1 : j*n_j;
        dataB(idx) = mu(type) + sig(type) * randn(n_j, 1);
        idsB(idx) = j;
    end

    % Visualization
    figure('Color', 'w', 'Position', [100 100 1000 600]);
    target_sub = 1;
    subDataA = dataA(idsA == target_sub);
    subDataB = dataB(idsB == target_sub);

    % Population plots
    subplot(2,2,1); histogram(dataA, -3:0.25:9, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    title('Population A: Shifters'); xlim([-3 8]); grid on;
    
    subplot(2,2,2); histogram(dataB, -3:0.25:9, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    title('Population B: Loyalists'); xlim([-3 8]); grid on;

    % Individual plots
    subplot(2,2,3); histogram(subDataA, -3:0.5:9, 'FaceColor', [0.8 0.2 0.2], 'EdgeColor', 'k');
    title(['Subject ' num2str(target_sub) ' (Shifter)']); xlim([-3 8]); grid on;
    
    subplot(2,2,4); histogram(subDataB, -3:0.5:9, 'FaceColor', [0.2 0.2 0.8], 'EdgeColor', 'k');
    title(['Subject ' num2str(target_sub) ' (Loyalist)']); xlim([-3 8]); grid on;
end

function results = cross_fit_demonstration(dataA, idsA, dataB, idsB, M)
    x0 = [0.5, 0, 5, 1, 1];
    lb = [0.001, -5, 0, 0.1, 0.1]; ub = [0.999, 5, 10, 5, 5];
    options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');

    % Fit to Shifter Data
    [resAA, fA] = fmincon(@(t) calc_negLL_A(t, dataA), x0, [], [], [], [], lb, ub, [], options);
    [resAB, fB] = fmincon(@(t) calc_negLL_B(t, dataA, idsA, M), x0, [], [], [], [], lb, ub, [], options);
    
    % Fit to Loyalist Data
    [resBA, fC] = fmincon(@(t) calc_negLL_A(t, dataB), x0, [], [], [], [], lb, ub, [], options);
    [resBB, fD] = fmincon(@(t) calc_negLL_B(t, dataB, idsB, M), x0, [], [], [], [], lb, ub, [], options);

    results.resAB = resAB; results.resBB = resBB; % Save for identity extraction
    
    % Console Output
    fprintf('\nSHIFTER DATA: Model A BIC %.1f vs Model B BIC %.1f\n', 2*fA + 5*log(length(dataA)), 2*fB + 5*log(length(dataA)));
    fprintf('LOYALIST DATA: Model A BIC %.1f vs Model B BIC %.1f\n', 2*fC + 5*log(length(dataB)), 2*fD + 5*log(length(dataB)));
end

function plot_identity_distributions(results, dataA, idsA, dataB, idsB, M)
    % Extract identities using Model B (The Identity Decoder)
    W_on_Shifters = get_identities(results.resAB, dataA, idsA, M);
    W_on_Loyalists = get_identities(results.resBB, dataB, idsB, M);

    figure('Color', 'w', 'Position', [150 150 1000 400]);
    
    subplot(1,2,1);
    histogram(W_on_Shifters, 0:0.05:1, 'FaceColor', 'r');
    title('Identities of Shifters (Model B perspective)');
    xlabel('Probability of being Type 1'); ylabel('Num Subjects');
    
    subplot(1,2,2);
    histogram(W_on_Loyalists, 0:0.05:1, 'FaceColor', 'b');
    title('Identities of Loyalists (Model B perspective)');
    xlabel('Probability of being Type 1');
end

function W = get_identities(theta, data, ids, M)
    p=theta(1); m1=theta(2); m2=theta(3); s1=theta(4); s2=theta(5);
    W = zeros(M, 1);
    for j = 1:M
        x = data(ids == j);
        l1 = log(p) + sum(log_normpdf(x, m1, s1));
        l2 = log(1-p) + sum(log_normpdf(x, m2, s2));
        W(j) = 1 / (1 + exp(l2 - l1));
    end
end

% --- Numerical Engines ---
function nll = calc_negLL_A(theta, data)
    lp1 = log(theta(1)) + log_normpdf(data, theta(2), theta(4));
    lp2 = log(1-theta(1)) + log_normpdf(data, theta(3), theta(5));
    nll = -sum(max(lp1,lp2) + log(exp(lp1-max(lp1,lp2)) + exp(lp2-max(lp1,lp2))));
end

function nll = calc_negLL_B(theta, data, ids, M)
    p=theta(1); m1=theta(2); m2=theta(3); s1=theta(4); s2=theta(5);
    total_logL = 0;
    for j = 1:M
        x = data(ids == j);
        l1 = log(p) + sum(log_normpdf(x, m1, s1));
        l2 = log(1-p) + sum(log_normpdf(x, m2, s2));
        total_logL = total_logL + max(l1, l2) + log(exp(l1-max(l1,l2)) + exp(l2-max(l1,l2)));
    end
    nll = -total_logL;
end

function v = log_normpdf(x, m, s)
    v = -0.5*((x-m)./s).^2 - log(s) - 0.5*log(2*pi);
end