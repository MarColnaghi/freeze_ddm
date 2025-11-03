function [est_kde] = extract_kde(bouts, points, export)

path = '/Users/marcocolnaghi/PhD/freeze_ddm/model_results/fitting_freezes/bsl/kde_spontaneous';
col = cmapper();
fh = figure('color','w','Position',[100, 100, 600, 400]);

x = 0:1/300:120;
[ f, xk ] = ksdensity(bouts.durations_s, x, 'Function','pdf', 'Bandwidth', 1/120);
[ F, ~  ] = ksdensity(bouts.durations_s, x, 'Function','cdf', 'Bandwidth', 1/120);

t0 = points.truncation;
F_t0 = interp1(xk, F, t0, 'linear', 'extrap');      % CDF at truncation
Z = max(1 - F_t0, eps);                              % normalizer (avoid /0)

fkde = zeros(size(f));
Fkde = zeros(size(F));

mask = xk >= t0;
fkde(mask) = f(mask) ./ Z;
Fkde(mask) = max( (F(mask) - F_t0) ./ Z, 0 );        % numerical safety

est_kde.xkde = xk;
est_kde.fkde = fkde;
est_kde.Fkde = Fkde;

% [est_kde.fkde, est_kde.xkde] = ksdensity(bouts.durations_s, 0:1/300:60, 'Function', 'pdf', 'Support','positive', 'Bandwidth', 1/10);
% [est_kde.Fkde, ~] = ksdensity(bouts.durations_s, 0:1/300:60, 'Function', 'cdf', 'Support','positive', 'Bandwidth', 1/10);
% est_kde.Fkde(est_kde.xkde < points.truncation) = 0;
% est_kde.fkde(est_kde.xkde < points.truncation) = 0;

histogram(bouts.durations_s, -1/120:1/60:300, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceColor', col.bsl)
hold on
plot(est_kde.xkde, est_kde.fkde, 'r--')
ax = gca;
apply_generic(ax)
xlim([-0.025 1.025]); xticks([0 1]);
ylim([-0.150 15.150])
xlabel('Immobility Duration (s)')
ax.YAxis.Visible = 'off';

if export
    exporter(fh, paths, 'kde_estimation.pdf')
end