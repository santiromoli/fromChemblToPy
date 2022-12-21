clear all, clc

%% build the regressors
                                                             #      1           2      3      4           5              6
load Xt_FGFR11.txt        %Xtraining --> ligands properties: ['QED Weighted', 'HBA', 'HBD', 'PSA', 'CX Acidic pKa', 'CX Basic pKa']
load yt_FGFR11.txt        %Ytraining --> IC50

xtotal = [Xt_FGFR11, yt_FGFR11];
xtotal = xtotal(all(~isnan(xtotal),2),:);

Xt = [xtotal(:,2) xtotal(:,4) xtotal(:,5) xtotal(:,6)];  #the most significant components

yt = [xtotal(:,7)];
yt = yt.^(1/14); % work in the power domain

startup      %call the GPML toolbox

SE = {@covSEard}; ell = 100; sf = 5*std(yt); hypSE = [log(ell)*ones(4,1); log(sf)];
CONS = {@covConst}; sf = 5*std(yt); hypCONS = log(sf);
cf = {'covSum',{SE,CONS}}; hyp.cov = [hypSE;hypCONS];
mf = {@meanSum, {@meanConst, @meanLinear}}; c = 0.0; hyp.mean = [c; zeros(4,1)];
lf = @likGauss; sn = std(yt); hyp.lik = log(sn);

% train the GPML (likGauss, meanZero)
[X, fX, i] = minimize(hyp, @gp, -100, @infExact, mf, cf, lf, Xt, yt);

%Infer the IC50 value of CHEMBL4083151-FGFR1
load Xtest.txt                  %CHEMBL4083151's properties

[m s2] = gp(X, @infExact, mf, cf, lf, Xt, yt, Xtest); #true IC50 value = 1.3

ii = 1;
ytest(ii) = m^(14)
hyps(ii,:) = [X.cov];
re = abs(1.3/ytest - 1);

##Loop

hyp.cov = [X.cov];
[X, fX, i] = minimize(hyp, @gp, -100, @infExact, mf, cf, lf, Xt, yt);
likelihood = X.lik;
hyps(ii+1,:) = [X.cov];
[m s2] = gp(X, @infExact, mf, cf, lf, Xt, yt, Xtest); #1.3
ytest(ii+1) = m^(14)

while (abs(ytest(ii)-ytest(ii+1)) > .01)
  hyp.cov = hyps(ii+1,:);
  [X, fX, i] = minimize(hyp, @gp, -100, @infExact, mf, cf, lf, Xt, yt);
  likelihood = X.lik;
  [m s2] = gp(X, @infExact, mf, cf, lf, Xt, yt, Xtest); #1.3
  ii = ii+1;
  ytest(ii+1) = m^(14)
  hyps(ii+1,:) = [X.cov];
endwhile
