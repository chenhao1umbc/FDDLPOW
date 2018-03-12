close all
clear
clc;
tic

addpath(genpath('./fddlow'))
addpath(genpath('./data'))
% load data
[Database]=load_data;

%% training dictionary
% for table 1 algorithm
K=100;
lbmd=1.5;
mu=0.1;
nu=1e3;
Q=16;% this is wq without negative
beta = 1;
SNR = 2000;
[opts]=loadoptions(K,lbmd,mu,Q,nu,beta, SNR);

% for table 2 algorithm
Dict_mix = FDDLOW_table2(Database.tr_data,Database.tr_label,opts);
save(opts.mixnm,'FDDLOW_mix','opts')

% this part will give the logistic regression result
trainLR_dif(opts, Database)


