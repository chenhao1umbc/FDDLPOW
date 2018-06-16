% table 2

close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))

% do traing or do crossvalidation
do_training = 1;
do_cv = 0;


% load data
mixture_n = 3; % mixture_n classes mixture, = 1,2,3
SNR = 2000;
pctrl.db = 10; % dynamic ratio is 0 3, 6, 10, 20 db
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
nu= 0.005 ;
beta = -1;
% nu = [ 0.003 0.005 0.007 0.01 0.03 0.05 0.07 0.1];
   
for ind1 = 1: length(nu)
    [opts]=loadoptions(K,lbmd,mu,Q,nu(ind1),beta, SNR);
    % for table 1 algorithm
    if do_training ==1
        Dict_mix = FDDLOW_table2(Database.tr_data,Database.tr_label,opts);
        if Dict_mix.iter > 30
            save(opts.mixnm,'Dict_mix','opts')
        end
    end
end

%% testing part
if do_cv ==1      
    result_nu = zeros(length(K),length(nu));
    sparsity_nu = zeros(length(K),length(nu));
    tr_sparsity_nu = zeros(length(K),length(nu));
    result_nuWEEK = zeros(length(K),length(nu));
    for ind1 = 1: length(nu)
        [opts]=loadoptions(K,lbmd,mu,Q,nu(ind1),beta, SNR);
        if exist(opts.mixnm, 'file')
            load(opts.mixnm)                            
            % run prep_ZF 
            if exist('Dict')==1
                Dict_mix = Dict;
            end
            Z = sparsecoding_mix_test(Dict_mix, Database, opts); %%%%% cv or test **************
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
            N_t = size(Database.test_mixlabel, 2); %%%%% cv or test **************
            opts.ln_test = N_t/featln;
            opts.equal = pctrl.equal;
            [acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts);

            result_nu(ind1) = acc_all
            result_nuWEEK(ind1) = acc_weak_av
            sparsity_nu(ind1) = mean(sum(Z ~= 0))/K;
            tr_sparsity_nu(ind1) = mean(sum(Dict_mix.Z ~= 0))/K;

            [SW,SB]=calcfisher(Dict_mix.Z,Database.tr_label,opts);
            fWZ=trace(W'*SW*W)-trace(W'*SB*W)+norm(W'*Dict_mix.Z,'fro')^2;              
        end 
    end
end

save('tb2_results','result_nu','result_nuWEEK','sparsity_nu','tr_sparsity_nu')
toc

