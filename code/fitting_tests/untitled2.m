idx = 4
figure; hold on;
plot(mu(idx,:))
yline(theta(idx))
xline(round(ts(idx) * 60))
plot(cum_surv(idx,:))