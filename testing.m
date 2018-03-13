clear 
clc

addpath(genpath('./fddlow'))
addpath(genpath('./data'))

% load('DDLMDmix4_k100_lmbd1.5_mu0.1_Q16_nu1000iter_100.mat')
load('Delta_FDDLOW_mix_k100_lmbd1.5_mu0.1_Q16_nu1000_beta0.1.mat')
% load('FDDLOW_mix_k100_lmbd1.5_mu0.1_Q16_nu1000_beta10.mat')
mixture_n = 2; % mixture_n classes mixture
[Database]=load_data(mixture_n);% the equal power mixture, 400 samples per combination
if exist('Dict')==1
    Dict_mix = Dict;
end 
Z = sparsecoding_mix_test(Dict_mix, Database, opts);

load('B_X_Y.mat')
% load('SNR2000B_X_Y_mix.mat')
[acc, acc_av] = lr_test(Dict_mix, Database, Z, B);