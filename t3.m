% table 3

close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))


% do traing or do crossvalidation
do_training = 0;
do_result = 1;
cv = 0; % validation or testing

for uuu = [0 3 6 10 20]
% load data
mixture_n = 3; % mixture_n classes mixture, = 1,2,3
SNR = 2000;
pctrl.db = uuu; % dynamic ratio is 0 3, 6, 10, 20 db
if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end

% the equal power mixture, 400 samples per combination
[Database]=load_data_new(mixture_n, SNR, pctrl);
% [Database]=load_data;
% Database = load_data_spectr(1);

%% training dictionary
% load settings
K = 100;
lbmd = 1e-4;
mu=1e-3;
Q=16;% this is wq without negative
SNR = 2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from table one we know that there are combinationes accuracy is above 0.99
% one is K = 100, lambda = 1e-4, mu = 1e-3
% another is K = 100, lambda = 1e-3, mu = 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu= 0.01 ;
beta = 100; % or 0.14
% beta = [1e-6 1e-5 1e-4 1e-3 1e-2 0.1 1 10 100 ];
   
for ind1 = 1: length(beta)
    [opts]=loadoptions(K,lbmd,mu,Q,nu,beta(ind1), SNR);
    % for table 1 algorithm
    if do_training ==1
        Dict_mix = FDDLOW_table3(Database.tr_data,Database.tr_label,opts);
        if Dict_mix.iter > 30
            save(opts.mixnm,'Dict_mix','opts')
        end
    end
end

%% testing part
if do_result ==1      
    result_beta = zeros(length(K),length(beta));
    sparsity_beta = zeros(length(K),length(beta));
    tr_sparsity_beta = zeros(length(K),length(beta));
    result_betaWEEK = zeros(length(K),length(beta));
    for ind1 = 1: length(beta)
        [opts]=loadoptions(K,lbmd,mu,Q,nu,beta(ind1), SNR);
        if exist(opts.mixnm, 'file')
            load(opts.mixnm)                            
            % run prep_ZF 
            if exist('Dict')==1
                Dict_mix = Dict;
            end
            if cv == 1;
                Z = sparsecoding_mix_cv(Dict_mix, Database, opts); %%%%% cv or test **************
            else
                Z = sparsecoding_mix_test(Dict_mix, Database, opts);
            end
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
            if cv == 1
                N_t = size(Database.cv_mixlabel, 2); %%%%% cv or test **************
            else 
                N_t = size(Database.test_mixlabel, 2); %%%%% cv or test **************
            end
            opts.ln_test = N_t/featln;
            opts.equal = pctrl.equal;
            [acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts);

            result_beta(ind1) = acc_all
            result_betaWEEK(ind1) = acc_weak_av
            sparsity_beta(ind1) = mean(sum(Z ~= 0))/K;
            tr_sparsity_beta(ind1) = mean(sum(Dict_mix.Z ~= 0))/K;

            [SW,SB]=calcfisher(Dict_mix.Z,Database.tr_label,opts);
            fWZ=trace(W'*SW*W)-trace(W'*SB*W)+norm(W'*Dict_mix.Z,'fro')^2;              
       
        end 
    end
end
end
 
save('tb3_results','result_beta','result_betaWEEK','sparsity_beta','tr_sparsity_beta')
toc