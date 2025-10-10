close all
clear vars
fitted_model = 'dddm2';
run = 'run03';

cd(fullfile('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/fitting_freezes/le', fitted_model, run))
model_fit = importdata('fit_results_dddm2.mat');
surrogate = importdata('surrogate.mat');

dx = .01;
xbin = (0:dx:dx*(ceil(10/dx)+10))';
xxi     = (0:10:20000)/1000;
h       = .05;
RTs{1,1} = surrogate.durations ./ 60;
RTs{2,1} = surrogate.durations ./ 60;

[RTD, TCM] = kreg_single(RTs, RTs, xxi, xbin, h, 0, 1000);

figure('color','w','Position',[100, 100, 900, 400]);
tiledlayout(1,3, 'TileSpacing', 'compact', 'Padding', 'compact')

nexttile
fill_between(xxi, zeros(1,length(RTD{1,2})), RTD{1,2} ,[],'FaceColor','#B9BFC2','FaceAlpha', 0.4,'LineStyle','none')
hold on
plot(xxi, RTD{1,2}, 'Color', '#B9BFC2', 'LineWidth',5)
