
idx_freeze = 2;
fh = figure('color','w','Position',[100, 100, 1300, 320]);
col = cmapper([], 1);
plot(mu(idx_freeze,:), 'Color', col.vars.sm(2,:), 'LineWidth', 2)
hold on
plot(1 - surv_probs(idx_freeze,:), 'Color', 'b', 'LineWidth', 2, 'DisplayName', 'P(non-break on frame)') % 
plot(cum_surv(idx_freeze,:), 'Color', 'g', 'LineWidth', 2, 'DisplayName', 'P(non-break until frame)') % Probability of non breaking up to that frame

xline(frames(idx_freeze), 'Label', 'Freeze Offset', 'FontSize', 16) % Fly broke from freezing 

scatter(frames(idx_freeze), cross_prob(idx_freeze), 80, 'o', 'MarkerFaceColor', 'b', 'DisplayName', 'P(non-break on frame)')
scatter(frames(idx_freeze) - 1, survival_prob(idx_freeze), 80, 'o', 'MarkerFaceColor', 'g', 'DisplayName', 'P(non-break until frame)')
legend
legend('Box', 'off', 'Location', 'southeastoutside')
apply_generic(gca)

plot(pdf, 'k')
sum(pdf)
scatter(frames(idx_freeze), survival_prob(idx_freeze) * cross_prob(idx_freeze), 80, 'o', 'MarkerFaceColor', 'k', 'DisplayName', 'P(non-break until frame)')
