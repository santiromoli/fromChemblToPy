% run startup.m prior
% activity estimator for TAS-120 against erbB-2
clear all, close all
pkg load statistics
rand('seed',0)

%% build the regressors

load Xt.txt        % Xtraining --> ligands properties: QED, CX Acidic pKa, CX Basic pKa,	HBA,	HBD,	PSA
load yt.txt        % Ytraining --> IC50
load Xtest.txt     % TAS-120's properties

% replace NaN values by mean
xm = nanmean(Xt);
for ii=1:6
  Xt(isnan(Xt(:,ii)),ii) = xm(ii);
  Xtest(isnan(Xtest(:,ii)),ii) = xm(ii);
end

% yt = log(yt); % work in the log domain
yt = yt.^(1/7); % work in the power domain

% subsample training set
ii = randperm(size(Xt,1));
n = 1000;
Xt = Xt(ii(1:n),:);
yt = yt(ii(1:n));

%retardo = 1; % delay

%xd = delay(xent,retardo);           % if desired
%xd2 = delay(xd,retardo);

%SE = {@covSEard};
%L = rand(6,1); sf = 2; hypSE = log([L;sf]);
%LIN = {'covLINard'};
%L = rand(6,1); hypLIN = (L);
%CONS = {@covConst};
%sf = 2;  hypCONS = log(sf);
%cf = {'covSum',{SE,LIN,CONS}};      %definition of covariance function
%hyp.cov = [hypSE;hypLIN;hypCONS];

cf = @covSEard; ell = 100; sf = 5*std(yt); hyp.cov = [log(ell)*ones(6,1); log(sf)];
mf = {@meanSum, {@meanConst, @meanLinear}}; c = 0.0; hyp.mean = [c; zeros(6,1)];
lf = @likGauss; sn = std(yt); hyp.lik = log(sn);

% train the GPML (likGauss, meanZero)
[X, fX, i] = minimize(hyp, @gp, -100, @infExact, mf, cf, lf, Xt, yt);
likelihood = X.lik % has to be <0
hyps = [X.cov];

%xd = delay(xval,retardo);
%xd2 = delay(xd,retardo);

[m s2] = gp(X, @infExact, mf, cf, lf, Xt, yt, Xtest);
m
