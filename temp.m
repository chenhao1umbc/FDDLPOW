clear 
clc

addpath(genpath('./fddlow'))
addpath(genpath('./data'))
table_weak = zeros(4,5);
table_av = zeros(4,5);
db_table = [3, 6, 10, 20, 40];

for ind1 = 1:4
    for ind2 = 1:5
        if ind1 ==1 || ind1 ==2
            load('DDLMDmix4_k100_lmbd1.5_mu0.1_Q16_nu1000iter_100.mat')%FDDLO
        else
            load('FDDLOW_mix_k100_lmbd1.5_mu0.1_Q16_nu1000_beta0.5.mat')
        end
        if ind1 ==1
            load('B_X_Y.mat')%FDDLO
        end
        if ind1 ==2 
            load('B_X_Y_pure.mat')%FDDLO
        end       
        if ind1 ==3
            load('SNR2000_beta0.5B_X_Y.mat')
        end
        if ind1 ==4
            load('beta0.5B_X_Y_pure.mat')
        end

        mixture_n = 3; % mixture_n classes mixture
        SNR = 2000;
        pctrl.equal = 0; % 1 means eqaul power, 0 non-equal
        pctrl.db = db_table(ind2); % dynamic ratio is 3, 6, 10, 20, 40db
        % the equal power mixture, 400 samples per combination
        [Database]=load_data(mixture_n, SNR, pctrl);
        if exist('Dict')==1
            Dict_mix = Dict;
        end 
        Z = sparsecoding_mix_test(Dict_mix, Database, opts);

        sparsity = mean(sum(Z ~= 0))
        [acc, acc_weak_av, acc_av] = lr_test(Dict_mix, Database, Z, B, pctrl);
        table_weak(ind1, ind2) = acc_weak_av;
        table_av(ind1, ind2) =acc_av;
    end
end