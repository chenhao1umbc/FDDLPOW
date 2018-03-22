% temp file to tune beta with equalizer

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
%             load('B_X_Y.mat')%FDDLO
        else
            load(['FDDLOW_mix_k100_lmbd1.5_mu0.1_Q16_nu1000_beta', num2str(t(ii)), '.mat'])
%             load(['SNR2000_beta', num2str(t(ii)), 'B_X_Y.mat'])
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
        
        W = Dict_mix.W;
        C = max(Database.tr_label);
        N = size(Database.tr_label,2);
        Nc = N / C;
        H3 = kron(eye(C),ones(Nc, 1)/Nc); % M = Z*H3
        M = Dict_mix.Z*H3;
        % zero forcing
        H = W'*M;
        result = pinv(H)*W'*aoos(Z,4,size(Z, 2));
        [~, labels_pre] = sort(result, 1, 'descend');
        
        opts.C = C; % 6 classes
        featln = Database.featln;
        opts.n = Database.N_c;
        opts.Ncombs = max(Database.cv_mixlabel);
        N_t = size(Database.test_mixlabel, 2); % test signal length
        opts.ln_test = N_t/featln;

        [acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts)

        record_weak_av(ii, jj) = acc_weak_av; 
        record_all(ii, jj) = acc_all; 
    end
end