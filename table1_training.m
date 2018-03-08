close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
% load data
[Database]=load_data;

%% training dictionary
% for Table 1 algorithm
K=100;
lambda = 1; %lambda=10;
mu = 0.1; %mu=1;
Q = 16;
[opts]=loadoptions(K,lambda,mu,Q);

% for table 1 algorithm
Dict = FDDLOW_table1(Database.tr_data, Database.tr_label, opts);
save(['.././', opts.mixnm],'Dict','opts')

