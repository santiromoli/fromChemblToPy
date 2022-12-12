#Activity estimator for FGFR-1 - CHEMBL4083151

clear all, clc

#build the regressors

load Xt_FGFR11.txt        #Xtraining --> ligands properties: QED, CX Acidic pKa,	CX Basic pKa,	HBA,	HBD,	PSA
load yt_FGFR11.txt        #Ytraining --> IC50

#remove all NaN inputs
xtotal = [Xt_FGFR11, yt_FGFR11];
xtotal = xtotal(all(~isnan(xtotal),2),:);

Xt = [xtotal(:,1) xtotal(:,2) xtotal(:,3) xtotal(:,4) xtotal(:,5) xtotal(:,6) xtotal(:,1).*xtotal(:,6)];
yt = [xtotal(:,7)];
#we work in the power domain in order to reduce squared distance of yt
yt = yt.^(1/14); 

startup      #call the GPML toolbox

cf = @covSEard; ell = 100; sf = 5*std(yt); hyp.cov = [log(ell)*ones(7,1); log(sf)];
mf = {@meanSum, {@meanConst, @meanLinear}}; c = 0.0; hyp.mean = [c; zeros(7,1)];
lf = @likGauss; sn = std(yt); hyp.lik = log(sn);

#train the GPML (likGauss, meanZero)
[X, fX, i] = minimize(hyp, @gp, -100, @infExact, mf, cf, lf, Xt, yt);
likelihood = X.lik % -1.5151 | -1.5165
hyps = [X.cov];

#infer the IC50 value of CHEMBL4083151-FGFR1

load Xtest.txt                  #CHEMBL4083151's properties

[m s2] = gp(X, @infExact, mf, cf, lf, Xt, yt, Xtest); #1.3

m^(14)
