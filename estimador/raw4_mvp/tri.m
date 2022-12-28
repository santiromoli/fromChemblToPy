clear all, clc

#build the regressors                     
load Xt_todo.txt        #Xtraining --> ligands properties: ['QED Weighted', 'HBA', 'HBD', 'PSA', 'CX Acidic pKa', 'CX Basic pKa' 
                                      #'cx_logd', 'cx_logp', 'alogp', 'mw_freebase', 'full_mwt', 'num_ro5_violations', 'heavy_atoms']
load yt_todo.txt        #Ytraining --> IC50

xtotal = [Xt_todo, yt_todo];
xtotal = xtotal(all(~isnan(xtotal),2),:); #take off all the NaN

Xt = [xtotal(:,7) xtotal(:,8) xtotal(:,9) xtotal(:,10) xtotal(:,11) xtotal(:,12) xtotal(:,13)];
Xt = (abs(Xt)).^(1/14);
yt = [xtotal(:,14)];
yt = yt.^(1/14); % work in the power domain

startup      %call the GPML toolbox
#load hypse.txt
se = {@covSEard}; ell = rand(7,1); #sf = 10*std(yt); hypse = [log(ell); log(sf)];
cons = {@covConst}; sf = 5*std(yt); #hypcons = log(sf);
load hypsecons.txt
hyp.cov = hypsecons;
cf = {'covSum',{se,cons}}; #hyp.cov = [hypse;hypcons];
mf = {@meanSum, {@meanConst, @meanLinear}}; c = 0.0; hyp.mean = [c; zeros(7,1)];
lf = @likGauss; sn = std(yt); hyp.lik = log(sn);

#[nlz dnlz] = minimize(hyp, @gp, -100, @infExact, mf, cf, lf, Xt, yt)

#train the GPML (likGauss, meanZero)
[X, fX, i] = minimize(hyp, @gp, -100, @infExact, mf, cf, lf, Xt, yt);

#Infer the IC50 value of CHEMBL4083151-FGFR1
load Xtest.txt                  #CHEMBL4083151's properties
Xtest = (Xtest).^(1/14);
[m s2] = gp(X, @infExact, mf, cf, lf, Xt, yt, Xtest); #1.3

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

