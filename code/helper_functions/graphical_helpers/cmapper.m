function col = cmapper(~, quantiles)
% CMAPPER  Centralized color manager
%   col = cmapper([], quantiles)
%
%   quantiles controls the number of discrete levels (+1 for background)

if nargin < 2 || isempty(quantiles)
    quantiles = 4;
end

nq = quantiles + 1;

%% =====================================================================
%% Base palettes
%% =====================================================================

pal.purple  = cbrewer2('Purples', nq);
pal.blue    = cbrewer2('Blues',   nq);
pal.red     = cbrewer2('Reds',    nq);
pal.orange  = cbrewer2('Oranges', nq);
pal.green   = cbrewer2('Greens',  nq);
pal.grey    = cbrewer2('Greys',   nq);

%% =====================================================================
%% Empirical / generic
%% =====================================================================

col.empirical.col = [cbrewer2('Set2',3); cbrewer2('Dark2',3)];
col.empirical.col = col.empirical.col([1,4,2,5,3,6],:);
col.empirical.gray = [[240, 167, 164]; [226, 155, 152]; [213, 133, 130];[240, 167, 164]; [226, 155, 152]; [213, 133, 130]] /255;
col.empirical.gray = col.empirical.gray([1,4,2,5,3,6],:);
col.empirical.gray2 = [[.8, .8, .8]; [.75, .75, .75]; [.7, .7, .7];[.8, .8, .8]; [.75, .75, .75]; [.7, .7, .7]];
col.empirical.gray2 = col.empirical.gray2([1,4,2,5,3,6],:);
% col.empirical.reds = [pal.red([2 3 4], :); pal.red([2 3 4], :) - 0.1];
% col.empirical.reds = col.empirical.reds([1,4,2,5,3,6],:);

col.pmix_map = 'L17';

%% =====================================================================
%% Process / model colors (single hex values)
%% =====================================================================

col.single_proc = '#E69665';
col.double_proc = '#CE55C2';
col.triple_proc = '#F1A208';

col.processes.short   = '#E36588';
col.processes.long    = '#648DE5';
col.processes.contam = '#E69665';

col.period.bsl  = col.processes.contam;
col.period.loom = '#778EE4';

%%%

col.param.theta = '#0000FF';
col.param.mu = '#DF3365';

%% =====================================================================
%% Model colors
%% =====================================================================

col.model.e0   = '#ED7B84';
col.model.w0   = '#E2DE6F';
col.model.g0   = '#7EB77F';
col.model.p0   = '#749FB6';

col.model.de1   = '#7DEE9F';
col.model.dde1  = '#3AD067';
col.model.dkde1 = '#83C0E9';
col.model.ddkde1 = '#1690E2';

col.closedform = '#AF8BD1';
col.extremadetection = '#DF2CC7';

%% =====================================================================
%% Variable → palette mapping
%% =====================================================================
% RULE:
%   - If a variable is listed here → use its palette
%   - Otherwise → use col.vars.default (purple)

col.vars.default = pal.purple;   % <---- NEW GENERIC FALLBACK 🔥

% Explicit mappings
col.vars.ln   = pal.purple;
col.vars.nloom = pal.purple;
col.vars.nloom1 = pal.purple;

col.vars.ls   = pal.blue;
col.vars.sloom = pal.blue;
col.vars.sloom1 = pal.blue;

col.vars.sm   = pal.red;
col.vars.avg_sm = pal.red;
col.vars.avg_sm_freeze_norm = pal.red;
col.vars.am_post = pal.red;
col.vars.sm_post = pal.red;

col.vars.sm_pre  = pal.orange;
col.vars.smp     = pal.orange;
col.vars.avg_ss  = pal.orange;
col.vars.time_since_last = pal.orange;
col.vars.onsets_loomaligned = pal.orange;
col.vars.onsets_loomaligned_norm = pal.orange;

col.vars.fs   = pal.green;
col.vars.sf   = pal.green;
col.vars.sf1  = pal.green;
col.vars.avg_fs_freeze_norm = pal.green;
col.vars.avg_fs_1s_norm = pal.green;

col.vars.moving_flies   = pal.grey;
col.vars.n_mov_flies1  = pal.grey;

col.vars.intercept   = [1 1 1];

%% =====================================================================
%% Visual aids / misc
%% =====================================================================

col.visual_aid.spontaneous = '#D664BE';
col.visual_aid.mu    = '#EE4266';
col.visual_aid.theta = '#007FFF';

col.opto = '#EE6464';

col.genotype = {'#784770', '#0F4D60', '#B34B4B'};

col.areas = cbrewer2('RdPu', 12);
col.areas(11,:) = [87 167 115]/255;

%% =====================================================================
%% Convenience colors
%% =====================================================================

col.loom = pal.blue(round(end/2),:);
col.bsl  = 'k';

%%

col.n_mov_flies = cbrewer2('RdBu', 5);

col.ll = cbrewer2('RdBu', 5);

%%
col.stationary_sm = '#FF595E';
col.timevarying_sm = '#1982C4';


%%
col.early_sm = '#FF715B';
col.late_sm = '#1EA896';


% Initialize the primary palette (limit Set2 to its max of 8)
col.pca = cbrewer2('Set2', min(nq, 8));

% Check if we need to append more colors
if nq > 8
    % Append colors from Dark2 for the remaining count
    extra_colors = cbrewer2('Dark2', nq - 8);
    col.pca = [col.pca; extra_colors];
end

col.pca = cbrewer2('Spectral', nq);
% col.pca = colorcet('I2','N', nq);

end

function cmap = get_var_cmap(col, varname)
%GET_VAR_CMAP  Safe access to variable colormap
if isfield(col.vars, varname)
    cmap = col.vars.(varname);
else
    cmap = col.vars.default;
end


end
