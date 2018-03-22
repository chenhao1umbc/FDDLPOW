clear 
clc

addpath(genpath('./fddlow'))
addpath(genpath('./data'))

t = [0, 0.1:0.2:0.7, 1, 1.3:0.2:1.7, 2:0.25:2.5, 2.7, 3, 5, 10];
tt = [3,6,10,20,40];
for ii = 1:length(t) % loop through beta
    for jj = 1:length(tt)
        if ii == 1
            load('DDLMDmix4_k100_lmbd1.5_mu0.1_Q16_nu1000iter_100.mat')%FDDLO
            load('B_X_Y.mat')%FDDLO
        else
            load(['FDDLOW_mix_k100_lmbd1.5_mu0.1_Q16_nu1000_beta', num2str(t(ii)), '.mat'])
            load(['SNR2000_beta', num2str(t(ii)), 'B_X_Y.mat'])
        end
        mixture_n = 2; % mixture_n classes mixture
        SNR = 2000;
        pctrl.equal = 0; % 1 means eqaul power, 0 non-equal
        pctrl.db = tt(jj); % dynamic ratio is 3, 6, 10, 20, 40db
        [Database]=load_data(mixture_n, SNR, pctrl);% the equal power mixture, 400 samples per combination
        if exist('Dict')==1
            Dict_mix = Dict;
        end 
        Z = sparsecoding_mix_test(Dict_mix, Database, opts);
        sparsity = mean(sum(Z ~= 0))
        [acc, acc_av] = lr_test(Dict_mix, Database, Z, B)

        record(ii, jj) = acc_av; 
    end
end