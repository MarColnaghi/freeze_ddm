
clearvars

model_st_single = importdata('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/fitting_freezes/le_new/ed5/run01_260325/model_results.mat');
model_tv_single = importdata('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/fitting_freezes/le_new/ed5/run02_260325/model_results.mat');
model_tv_double = importdata('/Users/marcocolnaghi/PhD/freeze_ddm/model_results/fitting_freezes/le_new/ded2/run04_260325/model_results.mat');

elbo.st_single = model_st_single.elbo;
elbo.tv_single = model_tv_single.elbo;
elbo.tv_double = model_tv_double.elbo;

figure()
bar(1:length(fieldnames(elbo)), [elbo.st_single(1), elbo.tv_single(1) elbo.tv_double(1)])