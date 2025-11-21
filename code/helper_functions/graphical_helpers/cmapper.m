function col = cmapper(string, quantiles)

if nargin < 2
    quantiles = [];
end

col.fits = [cbrewer2('Set2',3); cbrewer2('Dark2',3)];
col.fits = col.fits([1,4,2,5,3,6],:);
col.pmix_map = 'L17';
col.single_proc = '#E69665';
col.double_proc = '#CE55C2';
col.triple_proc =  '#F1A208';
col.gen_models.A = '#B8F4E2';
col.gen_models.B = '#F4C2CD';
col.gen_models.Aold = '#C1E5F5';
col.models.double = '#BFABCB';
col.models.singlecont = '#06A77D';
col.vars.intercept = [1 1 1];

col.vars.ln = cbrewer2('Purples', quantiles + 1);
col.vars.ls = cbrewer2('Blues', quantiles + 1);
col.vars.sm = cbrewer2('Reds', quantiles + 1);
col.vars.sm_post = cbrewer2('Reds', quantiles + 1);
col.vars.sm_pre = cbrewer2('Oranges', quantiles + 1);
col.vars.fs = cbrewer2('Greens', quantiles + 1);
col.vars.nloom1 = cbrewer2('Purples', quantiles + 1);
col.vars.sloom1 = cbrewer2('Blues', quantiles + 1);
col.vars.am_post1 = cbrewer2('Reds', quantiles + 1);
col.vars.am_pre1 = cbrewer2('Oranges', quantiles + 1);
col.vars.n_mov_flies1 = cbrewer2('Greys', quantiles + 1);
col.vars.sf1 = cbrewer2('Greens', quantiles + 1);
col.vars.nloom = cbrewer2('Purples', quantiles + 1);
col.vars.sloom = cbrewer2('Blues', quantiles + 1);
col.vars.am_post = cbrewer2('Reds', quantiles + 1);
col.vars.avg_sm = cbrewer2('Reds', quantiles + 1);
col.vars.sm = cbrewer2('Reds', quantiles + 1);
col.vars.smp = cbrewer2('Oranges', quantiles + 1);

col.vars.onsets_loomaligned = cbrewer2('Oranges', quantiles + 1);

col.vars.moving_flies = cbrewer2('Greys', quantiles + 1);
col.vars.sf = cbrewer2('Greens', quantiles + 1);
col.vars.avg_sm_freeze_norm = cbrewer2('Reds', quantiles + 1);
col.vars.avg_fs_freeze_norm = cbrewer2('Greens', quantiles + 1);
col.vars.onset_loomaligned = cbrewer2('Oranges', quantiles + 1);
col.vars.onsets_loomaligned_norm = cbrewer2('Oranges', quantiles + 1);

col.vars.avg_fs_1s_norm = cbrewer2('Greens', quantiles + 1);
col.visual_aid.spontaneous = '#D664BE';
col.visual_aid.mu = '#EE4266';
col.visual_aid.theta = '#007FFF';
col.opto = '#EE6464';
col.genotype = {'#784770', '#0F4D60', '#B34B4B'};
col.areas = cbrewer2('RdPu', 12); col.areas(11, :) = [87 167 115]/255;

col.processes.short = '#E36588';
col.processes.long = '#648DE5';
col.processes.contam = '#E69665';

col.loom = col.vars.sloom(round(end/2),:);

col.model.e0 = '#ED7B84';
col.model.e0first30frames = '#A51925';
col.model.e0first60frames = '#EB4250';
col.model.de1 = '#ED7B84';
col.model.de1_l30 = '#A51925';
col.model.de1_l60 = '#EB4250';
col.model.w0 = '#E2DE6F';
col.model.g0 = '#7EB77F';
col.model.p0 = '#749FB6';
col.model.dw0 = '#E2DE6F';
col.model.race0 = '#E23E6F';

col.model.de1 = '#7DEE9F';
col.model.dde1 = '#3AD067';
col.model.dkde1 = '#83C0E9';
col.model.ddkde1 = '#1690E2';

col.early_sm = '#FF715B';
col.late_sm = '#1EA896';

col.stationary_sm = '#FF595E';
col.timevarying_sm = '#1982C4';
col.extremadetection ='#DF2CC7';
col.closedform = '#AF8BD1';
col.bsl = 'k';

