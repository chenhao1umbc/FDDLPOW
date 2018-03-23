clear 
clc

addpath(genpath('./fddlow'))
addpath(genpath('./data'))

load('DDLMDmix4_k100_lmbd1.5_mu0.1_Q16_nu1000iter_100.mat')%FDDLO
load('B_X_Y.mat')%FDDLO

% load('FDDLOW_mix_k100_lmbd1.5_mu0.1_Q16_nu1000_beta3.mat')
% load('SNR2000_beta3B_X_Y.mat')

mixture_n = 2; % mixture_n classes mixture
SNR = 2000;
pctrl.equal = 1; % 1 means eqaul power, 0 non-equal
pctrl.db = 40; % dynamic ratio is 3, 6, 10, 20, 40db
[Database]=load_data(mixture_n, SNR, pctrl);% the equal power mixture, 400 samples per combination
if exist('Dict')==1
    Dict_mix = Dict;
end 
Z = sparsecoding_mix_test(Dict_mix, Database, opts);

sparsity = mean(sum(Z ~= 0))
[acc, acc_weak_av, acc_av] = lr_test(Dict_mix, Database, Z, B, pctrl)