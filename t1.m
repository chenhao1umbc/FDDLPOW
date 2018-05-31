% table 1 

close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))

% do traing or do crossvalidation
do_training = 0;
do_cv = 1;

% load data
mixture_n = 1; % mixture_n classes mixture, = 1,2,3
SNR = 2000;
pctrl.equal = 1; % 1 means eqaul power, 0 non-equal
pctrl.db = 3; % dynamic ratio is 3, 6, 10, 20, 40db

% the equal power mixture, 400 samples per combination
[Database]=load_data_new(mixture_n, SNR, pctrl);
% [Database]=load_data;
% Database = load_data_spectr(1);

%% training dictionary
% load settings
K = 100;
lbmd = 0.001;
mu=0.1;
nu= 1e3;
Q=16;% this is wq without negative
beta = 1;
SNR = 2000;

K = [100, 150, 200, 250];
lbmd = [0.001 0.005 0.01 1e-4];
mu = [1, 0.1, 0.01, 0.001 0.0001];

    
for ind1 = 1: length(K)
    for ind2 = 1: length(lbmd)
        for ind3= 1: length(mu)
            [opts]=loadoptions(K(ind1),lbmd(ind2),mu(ind3),Q,nu,beta, SNR);
            % for table 1 algorithm
            if do_training ==1
                Dict = FDDLOW_table1(Database.tr_data,Database.tr_label,opts);
                if Dict.iter > 80
                    save(opts.Dictnm,'Dict','opts')
                end
            end 

        end
    end
end

%% testing part
if do_cv ==1      
    result_K_lambda_mu = zeros(length(K),length(lbmd),length(mu));
    sparsity_K_lambda_mu = zeros(length(K),length(lbmd),length(mu));
    for ind1 = 1: length(K)
        for ind2 = 1: length(lbmd)
            for ind3= 1: length(mu)
                [opts]=loadoptions(K(ind1),lbmd(ind2),mu(ind3),Q,nu,beta, SNR);
                if exist(opts.Dictnm, 'file')
                    load(opts.Dictnm,'Dict','opts')
                end               
                % run prep_ZF 
                if exist('Dict')==1
                    Dict_mix = Dict;
                end
                Z = sparsecoding_mix_cv(Dict_mix, Database, opts);
                W = Dict_mix.W;
                C = max(Database.tr_label);
                N = size(Database.tr_label,2);
                Nc = N / C;
                opts.C = C; % 6 classes
                featln = Database.featln;
                opts.n = Database.N_c;                
                H3 = kron(eye(C),ones(Nc, 1)/Nc); % M = Z*H3
                M = Dict_mix.Z*H3;
                % zero forcing
                H = W'*M;
                result = pinv(H)*W'*aoos(Z,featln,size(Z, 2));
                [~, labels_pre] = sort(result, 1, 'descend');

                opts.Ncombs = max(Database.cv_mixlabel);
                N_t = size(Database.cv_mixlabel, 2); % test signal length
                opts.ln_test = N_t/featln;
                opts.equal = pctrl.equal;
                [acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts);
                
                result_K_lambda_mu(ind1, ind2, ind3) = acc_all;
                sparsity_K_lambda_mu(ind1, ind2, ind3) = mean(sum(Z ~= 0))/K(ind1);
                
                [SW,SB]=calcfisher(Dict_mix.Z,Database.tr_label,opts);
                fWZ=trace(W'*SW*W)-trace(W'*SB*W)+norm(W'*Dict_mix.Z,'fro')^2;
                
                mean(sum(Z ~= 0))/K(ind1)
                mean(sum(Dict_mix.Z ~= 0))/K(ind1)
                

                opts.lambda1*sum(abs(Dict_mix.Z(:)))
                opts.mu*fWZ
                
            end
        end
    end
end

toc


