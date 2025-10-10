function [RTD,TCM] = kreg_single(X,Y,xxi,xbin,h,RTmin,B)
% calculate RTDs and TDAs for each condition with bootstrap
% the number of conditions is inferred from the size of the cell RTs
% --- INPUTS ---
% - RTs: cell of RTs per condition (each new dimension is a new condition)
%       IMPORTANT: The first dimension MUST be binary, and it will be the
%       dimension for which the TCMometric is calculated (it could be
%       outcome, choice or other)
% - xxi: values of RT for which the TDAs and RTDs are calculated
% - h: kernel bandwith
% - RTmin: min RT (to account for boundary effects)
% - B: number of bootstrap samples
% for the smoothed bootstrap method, see:
% https://stats.stackexchange.com/questions/155488/simulate-from-kernel-density-estimate-empirical-pdf

%%

resamp = 1;
if B<=1
    B = 1;
    resamp = 0;
end

% h3  = h*1.2;

qs = [0.025,0.975]; %quantiles to estimate CIs
sq = @(x) squeeze(x);
Kh  = @(y) 3/4*(1-(y/h).^2)/h;
hk = linspace(0,h,500);
cf = zeros(size(hk));
for k=1:length(hk)
    cf(k) = 1./(.5 + integral(Kh,0,hk(k)));
end
cf_itp = griddedInterpolant(hk,cf,'spline','nearest');

epa = @(x) 3/4*(1-x.^2);
f_itp = griddedInterpolant(-1:.01:1,epa(-1:.01:1),'spline','nearest');
 
% y = (.1:dx:dx*(ceil(1/dx)+10))';
% Kn = cell(length(x),length(y));
z1   = (xxi-xbin)/h;
f1   = sparse(f_itp(z1)); %these are the individual kernels for each value of x
% z2   = (xxi-y)/h;
% f2   = sparse(f_itp(z2)'); %these are the individual kernels for each value of y
% for k1=1:length(x)
%     for k2=1:length(y)
%         Kn{k1,k2} = f1(k1,:).*f2(:,k2);
%     end
% end
% whos Kn    

n = size(X);
N = length(X(:));

temp = num2cell([B,length(xxi),n]);
Znum_BA = zeros([temp{:}]);
Zden_BA = zeros([temp{:}]);

for ki=1:N
    RTx     = X{ki};
    L = length(RTx);
    x_m = mean(RTx(:,1));
    x_s = std(RTx(:,1));
    Znum   = zeros(B,length(xxi));
    Zden   = zeros(B,length(xxi));
    X_i = RTx;%sortrows(RTx,1);
    Y_i = Y{ki};
    
    for k1=1:B
        [znum,zden] = custom_kreg(xbin,X_i,Y_i,h,RTmin,f1,f_itp,cf_itp,L,resamp,x_m,x_s);
        Znum(k1,:) = znum;
        Zden(k1,:) = zden;
    end
    
    temp = num2cell(n);
    [temp{:}] = ind2sub(n,ki);
    
    tempk = [{':'},{':'},temp];
    Znum_BA(tempk{:}) = Znum;
    Zden_BA(tempk{:}) = Zden;

        
end

ZTCM_BA = Znum_BA./Zden_BA;
ZCDF_BA = cumtrapz(xxi,Zden_BA,2);

if N==1
    TCM_q1 = sq(quantile(ZTCM_BA,qs(2),1))';
    TCM_q2 = sq(quantile(ZTCM_BA,qs(1),1))';
    TCM_md = sq(nanmedian(ZTCM_BA,1))';
    RTD_q1 = sq(quantile(Zden_BA,qs(2),1))';
    RTD_q2 = sq(quantile(Zden_BA,qs(1),1))';
    RTD_md = sq(nanmedian(Zden_BA,1))';
    CDF_q1 = sq(quantile(ZCDF_BA,qs(2),1))';
    CDF_q2 = sq(quantile(ZCDF_BA,qs(1),1))';
    CDF_md = sq(nanmedian(ZCDF_BA,1))';
else
    TCM_q1 = sq(quantile(ZTCM_BA,qs(2),1));
    TCM_q2 = sq(quantile(ZTCM_BA,qs(1),1));
    TCM_md = sq(nanmedian(ZTCM_BA,1));
    RTD_q1 = sq(quantile(Zden_BA,qs(2),1));
    RTD_q2 = sq(quantile(Zden_BA,qs(1),1));
    RTD_md = sq(nanmedian(Zden_BA,1));
    CDF_q1 = sq(quantile(ZCDF_BA,qs(2),1));
    CDF_q2 = sq(quantile(ZCDF_BA,qs(1),1));
    CDF_md = sq(nanmedian(ZCDF_BA,1));
end

TCM_up = abs(TCM_q1-TCM_md);
TCM_dn = abs(TCM_q2-TCM_md);

RTD_up = abs(RTD_q1-RTD_md);
RTD_dn = abs(RTD_q2-RTD_md);

CDF_up = abs(CDF_q1-CDF_md);
CDF_dn = abs(CDF_q2-CDF_md);

TCM = cell(3,N);
RTD = cell(3,N);
CDF = cell(3,N);
for k=1:N
    TCM{1,k} = TCM_md(:,k);
    TCM{2,k} = TCM_up(:,k);
    TCM{3,k} = TCM_dn(:,k);
    RTD{1,k} = RTD_md(:,k);
    RTD{2,k} = RTD_up(:,k);
    RTD{3,k} = RTD_dn(:,k);
    CDF{1,k} = CDF_md(:,k);
    CDF{2,k} = CDF_up(:,k);
    CDF{3,k} = CDF_dn(:,k);
end



end



function [RTD_num,RTD_den] = custom_kreg(x,RTs,Ys,h,RTmin,f,f_itp,cf_itp,L,resamp,x_m,x_s)
% calculate RTDs and TDA using kernels
% -t: points for which to calculate the KDE
% -X_i: sample
% -h: kernel bandwidth
% -RTmin
% -RTmax
% -f_itp: kernel interpolant
% -cf_itp: kernel normalization constant interpolant 
% -x_m: mean of sample
% -x_s: std of sample
% -L: length of sample

corvar = 0; %use variance correction
% https://stats.stackexchange.com/questions/155488/simulate-from-kernel-density-estimate-empirical-pdf

%% resample data
if resamp==1
    samp    = randsample(L,L,true);
    r1      = median(2*rand(L,3),2)-1;
    r2      = median(2*rand(L,3),2)-1;
    X_i     = RTs(samp,:);
    Y_i     = Ys(samp);
    X_i     = X_i + h*r1;
%     Y_i     = Y_i + h*r2;
    
%     if corvar==1
%         X_i(:,1) = x_m + (X_i(:,1) - x_m + h*r)./(sqrt(1+1/5*(h/x_s).^2));
%     else
%         X_i(:,1) = X_i(:,1) + h*r;
%     end
%     
%     % resample the RTs that came out negative
%     negX_lg = X_i(:,1)<=0; %logical vector of negative RTs
%     negX_in = find(negX_lg);
%     negL    = length(negX_in);
%     for kn=1:negL
%         while X_i(negX_in(kn),1)<=0
%             r  = median(2*rand(1,3),2)-1;
%             X_i(negX_in(kn),1) = RTs(samp(negX_in(kn)),1);
%             if corvar==1
%                 X_i(negX_in(kn),1) = x_m + (X_i(negX_in(kn),1) - x_m + h*r)./(sqrt(1+1/5*(h/x_s).^2));
%             else
%                 X_i(negX_in(kn),1) = X_i(negX_in(kn),1) + h*r;
%             end
%         end
%     end
else
    X_i = RTs;
    Y_i = Ys;
end

%% initialize and obtain kernels

% dx = .001;
% x = (0:dx:dx*(ceil(max(X_i(:,1))/dx)+10))'; %X_i will be binned for speedup
% x = (0:dx:dx*(ceil(.9/dx)+10))';

% T   = length(t);

%% to account for boundary effects

cf = ones(size(x));
dL = x-RTmin;

ind1 = (dL<0);
ind2 = (dL==0);
ind3 = (dL<h);

cf(ind1) = 0;
cf(ind2) = 2;
cf(ind3) = cf_itp(dL(ind3));

%% new method

[hc,~,bins] = histcounts(X_i,[x;x(end)+x(end)-x(end-1)]);
n = (1:length(hc))';
fun = @(ns) sum(Y_i(bins==ns));
hcy = arrayfun(fun,n).*cf;
hcx = hc'.*cf;


%% calculate KDEs

RTD_num = sum(hcy.*f)/(L*h);
RTD_den = sum(hcx.*f)/(L*h);


end

