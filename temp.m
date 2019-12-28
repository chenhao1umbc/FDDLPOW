% table 1 

close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././FDDLPOW'))
% do traing or do crossvalidation
do_training = 1;
do_cv = 0;

% load data
mixture_n = 1; % mixture_n classes mixture, = 1,2,3
pctrl.db = 0; % dynamic ratio is 0 3, 6, 10, 20 db
if pctrl.db == 0     pctrl.equal = 1; else    pctrl.equal = 0; end
cvortest = 1;  % 1 means cv, 0 means test

%% training dictionary
% load settings
K = 100;
lbmd = 1e-4;
mu=1e-3;
nu= 1e3;
beta = 1;
Q=8;% this is wq without negative
SNR = 2000;

% K = [1000];
% lbmd = [0.01];
% mu = [0.001, 0.0001];
% SNR = [2000, 20, 0, -5, -10, -20];
% Q = [16, 32, 48, 64, 80, 96];

[Database]=load_data_new(mixture_n, SNR, pctrl);

for ind1 = 1: length(K)
for ind2 = 1: length(lbmd)
for ind3= 1: length(mu)
    [opts]=loadoptions(K(ind1),lbmd(ind2),mu(ind3),Q,nu,beta, SNR)
    % for table 1 algorithm
    if do_training ==1
        Dict = FDDLOW_table1(Database.tr_data,Database.tr_label,opts);
        save(opts.Dictnm,'Dict','opts')
    end
end 
end
end