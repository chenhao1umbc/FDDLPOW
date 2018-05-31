clear 
clc

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))

% load('DDLMDmix4_k100_lmbd1.5_mu0.1_Q16_nu1000iter_100.mat') %FDDLO
% load(['rng',num2str(ii), 'FDDLOW_mix_k100_lmbd1.5_mu0.1_Q16_nu1000_beta1.7.mat'])
% load('FDDLOW_mix_k100_lmbd1.5_mu0.1_Q16_nu1000_beta1.7.mat')
% load('FDDLOW_mix_k100_lmbd1.5_mu0.1_Q32_nu1000_beta1.7.mat')

load('FDDLOW_mix_k100_lmbd1.5_mu0.1_Q16_nu1000_beta1.mat')

mixture_n = 2; % mixture_n classes mixture
SNR = 2000;
pctrl.equal = 1; % 1 means eqaul power, 0 non-equal
pctrl.db = 10; % dynamic ratio is 3, 6, 10, 20, 40db
% the equal power mixture, 400 samples per combination
[Database]=load_data_new(mixture_n, SNR, pctrl);

run prep_ZF 
[acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts)