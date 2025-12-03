% Load the table first. We will take advantage of an already existing
% dataset.
threshold_imm = 1; threshold_mob = 1; threshold_pc = 4; id_code = sprintf('imm%d_mob%d_pc%d', threshold_imm, threshold_mob, threshold_pc);
paths = path_generator('folder', 'fitting_freezes/bsl/kde_spontaneous', 'bouts_id', id_code);
load(fullfile(paths.dataset, 'bouts.mat'));
bouts_proc = data_parser_new(bouts, 'type', 'immobility', 'period', 'bsl', 'window', 'all');
points.censoring = 10.5;
points.truncation = min(bouts_proc.durations_s);

col = cmapper();
fh = figure('color','w','Position',[100, 100, 600, 400]);

x = 0:1/60:120;
[ f, xk ] = ksdensity(bouts_proc.durations_s, x, 'Function', 'pdf');
[ F, ~  ] = ksdensity(bouts_proc.durations_s, x, 'Function', 'cdf');

% t0 = points.truncation;
% F_t0 = interp1(xk, F, t0, 'linear', 'extrap');      % CDF at truncation
% Z = max(1 - F_t0, eps);                              % normalizer (avoid /0)
% 
% fkde = zeros(size(f));
% Fkde = zeros(size(F));
% 
% mask = xk >= t0;
% fkde(mask) = f(mask) ./ Z;
% Fkde(mask) = max( (F(mask) - F_t0) ./ Z, 0 );        % numerical safety

est_kde.xkde = xk;
est_kde.fkde = f;
est_kde.Fkde = F;
area_pdf = trapz(est_kde.xkde, est_kde.fkde);  % expect ~1

histogram(bouts_proc.durations_s, 1/120:1/60:300, 'Normalization', 'pdf', 'EdgeColor', 'none', 'FaceColor', col.bsl)
hold on
plot(est_kde.xkde, est_kde.fkde, 'r--')
ax = gca;
apply_generic(ax)
xlim([-0.025 1.025]); xticks([0 1]);
ylim([-0.150 15.150])
xlabel('Immobility Duration (s)')
ax.YAxis.Visible = 'off';

if saving
    mkdir(fullfile(paths.results, id_code))
    cd(fullfile(paths.results, id_code))
    save('kde_estimates_bsl.mat', 'est_kde')
end