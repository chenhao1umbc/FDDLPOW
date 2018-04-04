clear 
clc

addpath(genpath('./fddlow'))
addpath(genpath('./data'))
table_weak = zeros(2,5);
table_av = zeros(2,5);
db_table = [3, 6, 10, 20, 40];
for ind1 = 1:2
    for ind2 = 1:5
        if ind1  ==1
            load('DDLMDmix4_k100_lmbd1.5_mu0.1_Q16_nu1000iter_100.mat') %FDDLO
        else 
            load('FDDLOW_mix_k100_lmbd1.5_mu0.1_Q16_nu1000_beta1.7.mat')
        end

        mixture_n = 3; % mixture_n classes mixture
        SNR = 2000;
        pctrl.equal = 0; % 1 means eqaul power, 0 non-equal
        pctrl.db = db_table(ind2); % dynamic ratio is 3, 6, 10, 20, 40db
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
        opts.equal = pctrl.equal;

        [acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts)
        table_weak(ind1, ind2) = acc_weak_av;
        table_av(ind1, ind2) = acc_all;
    end
end