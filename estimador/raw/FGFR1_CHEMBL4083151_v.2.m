#Activity estimator for FGFR-1 - CHEMBL4083151

clear all, clc

#build the regressors

load Xt_FGFR11.txt        #Xtraining --> ligands properties: ['QED Weighted', 'HBA', 'HBD', 'PSA', 'CX Acidic pKa', 'CX Basic pKa']
load yt_FGFR11.txt        #Ytraining --> IC50

xtotal = [Xt_FGFR11, yt_FGFR11];
xtotal = xtotal(all(~isnan(xtotal),2),:);

Xt = [xtotal(:,1).*xtotal(:,4) xtotal(:,2) xtotal(:,3) xtotal(:,5) xtotal(:,6)];
yt = [xtotal(:,7)];
yt = yt.^(1/14);                                                                 % work in the power domain

startup      %call the GPML toolbox

#feval(@covSEard) e.g. to ask for the number of hyp
SE = {@covSEard}; ell = 100; sf = 5*std(yt); hypSE = [log(ell)*ones(5,1); log(sf)];
NO = {'covNoise'}; sf = 5*std(yt); hypNO = log(sf);                                    #.5*std?
LIN = {'covLINard'}; ell = 100; hypLIN = [log(ell)*ones(5,1)];
CONS = {@covConst}; sf = 5*std(yt); hypCONS = log(sf);
cf = {'covSum',{SE,NO,LIN,CONS}}; hyp.cov = [hypSE;hypNO;hypLIN;hypCONS];
mf = {@meanSum, {@meanConst, @meanLinear}}; c = 0.0; hyp.mean = [c; zeros(5,1)];
lf = @likGauss; sn = std(yt); hyp.lik = log(sn);

#train the GPML (likGauss, meanZero)
[X, fX, i] = minimize(hyp, @gp, -100, @infExact, mf, cf, lf, Xt, yt)
likelihood = X.lik
hyps = [X.cov];

%Infer the IC50 value of CHEMBL4083151-FGFR1

load Xtest.txt                  %CHEMBL4083151's properties

[m s2] = gp(X, @infExact, mf, cf, lf, Xt, yt, Xtest); #1.3

m^(14)

#hyp.cov = hyps;                   #re-initialize the optimization

#[X, fX, i] = minimize(hyp, @gp, -100, @infExact, mf, cf, lf, Xt, yt);
#likelihood = X.lik
#hyps = [X.cov];

#[m s2] = gp(X, @infExact, mf, cf, lf, Xt, yt, Xtest); #1.3
#m^(14)

