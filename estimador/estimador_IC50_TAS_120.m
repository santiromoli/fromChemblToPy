% run startup.m prior
% activity estimator for TAS-120 against erbB-2
clear all, close all
clc

%% build the regressors

load Xt.txt        %Xtraining --> ligands properties: QED, CX Acidic pKa,	CX Basic pKa,	HBA,	HBD,	PSA
load yt.txt        %Ytraining --> IC50
load Xtest.txt     %TAS-120's properties

%retardo = 1; %delay

%xd = delay(xent,retardo);           %if desired
%xd2 = delay(xd,retardo);

startup      %call the GPML toolbox
meanfunc = {@meanZero}; %this function has not hyp.mean
%feval(covfunc{:}) % how many hyperparametrs has meanfunc, covfunc or likfunc
SE = {@covSEard};
L = rand(6,1); sf = 2; hypSE = log([L;sf]);
LIN = {'covLINard'};
L = rand(6,1); hypLIN = (L);
CONS = {@covConst};
sf = 2;  hypCONS = log(sf); %hypCONS = 0;
NO = {'covNoise'};
sn = .2;  hypNO = log(sn); %hypNO = 0;
ovfunc = {'covSum',{SE,NO,LIN,CONS}};      %definition of covariance function
hyp.cov = [hypSE;hypNO;hypLIN;hypCONS];

%feval(likfunc{:})
likfunc = {@likGauss};       %choose the likelihood function
sn = 0.1;
hyp.lik = log(sn);

%Train the GPML
[X, fX, i] = minimize(hyp, @gp, -1500, @infExact, meanfunc, covfunc, likfunc, Xt, yt);
likelihood = X.lik %has to be <0
hyps = [X.cov];

%Infer the IC50 value of TAS-120

load Xtest.txt

%xd = delay(xval,retardo);
%xd2 = delay(xd,retardo);

[m s2] = gp(X, @infExact, meanfunc, covfunc, likfunc, Xt, yt, Xtest);

m
