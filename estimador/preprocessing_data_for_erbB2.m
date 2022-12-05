% Preprocessing data from 'retrieve_data_futibatinib.py'
% Test case: Futibatinib against erbB2

load XeYt.mat

yt = erbB2_IC50_value';
Xt = erbB2';

Xtest = [0.51	NA	2.58	8	1	108.39];

save yt.txt yt -ascii
save Xt.txt Xt -ascii
save Xtest.txt Xtest -ascii
