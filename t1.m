% table 1 

close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././FDDLPOW'))
% do traing or do crossvalidation
do_training = 0;
do_cv = 1;
mixture_n = 1; % mixture_n classes mixture, = 1,2,3
pctrl.db = 0; % dynamic ratio is 0 3, 6, 10, 20 db
pctrl.if2weak = 0; % if 2 weak components in mixture of 3
if pctrl.db == 0     pctrl.equal = 1; else    pctrl.equal = 0; end
cvortest = 1;  % 1 means cv, 0 means test

%% training dictionary
% load settings
K = 50;
lbmd = 0.01;
mu=0.001;
nu= 1e3;
beta = 1;
Q= 45;% this is wq without negative
SNR = 0;

% K = [50, 100, 150, 200, 250];
% lbmd = [0.1, 0.01, 0.001, 1e-4];
% mu = [1, 0.1, 0.01, 0.001 0.0001];
% SNR = [2000, 20, 0, -5, -10, -20];
% Q = [16, 32, 48, 64, 80, 96]; % [50, 45, 40, 30, 15, 8]

[Database]=load_data_new(mixture_n, SNR, pctrl);

for ind1 = 1: length(K)
for ind2 = 1: length(lbmd)
for ind3= 1: length(mu)   
    % for table 1 algorithm
    if do_training ==1
        [opts]=loadoptions(K(ind1),lbmd(ind2),mu(ind3),Q,nu,beta, SNR)
        Dict = FDDLOW_table1(Database.tr_data,Database.tr_label,opts);
        save(opts.Dictnm,'Dict','opts')
    end
end 
end
end

%% CV/testing part
for f = 1:5
[Database]=load_data_new(mixture_n, SNR, pctrl);
if do_cv ==1      
% result_K_lambda_mu = zeros(length(K),length(lbmd),length(mu));
% sparsity_K_lambda_mu = zeros(length(K),length(lbmd),length(mu));
% tr_sparsity_K_lambda_mu = zeros(length(K),length(lbmd),length(mu));
% result_K_lambda_muWEEK = zeros(length(K),length(lbmd),length(mu));

% for Q = [30, 32]
for ind1 = 1: length(K)
for ind2 = 1: length(lbmd)
for ind3 = 1: length(mu)
    [opts]=loadoptions(K(ind1),lbmd(ind2),mu(ind3),Q,nu,beta, SNR);
    if exist(opts.Dictnm, 'file') load(opts.Dictnm,'Dict','opts'), else break; end
    % run prep_ZF 
    if exist('Dict')==1    Dict_mix = Dict; opts,   end
    Z = sparsecoding(Dict, Database, opts, mixture_n, cvortest);
    Z = aoos(Z,Database.featln, size(Z, 2));
    
    % KNN classifier
    if Database.N_c == 1        
        Xtestorcv = Dict.W'*Z;
        Xtr = Dict.W'*Dict.Z;%aoos(Dict.Z,Database.featln, size(Dict.Z, 2));
        % KNN classifier
        acc_knn = myknn(Xtr, Xtestorcv, Database, cvortest) % k = 5 
    end

%     run calc_M
%     % zero forcing
%     H = W'*M;
%     r_zf = pinv(H)*W'*Z;
%     [~, labels_pre] = sort(r_zf, 1, 'descend');
%     
%     % matched filter
%     r_matched = H'*W'*Z;
%     [~, labels_pre_mf] = sort(r_zf, 1, 'descend');
%     
%     % calculate accuracy
%     [acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts);
%     [t, tt, ttt] = calc_labels(labels_pre_mf, opts);
% 
%     result_K_lambda_mu(ind1, ind2, ind3) = acc_all;
%     result_K_lambda_muWEEK(ind1, ind2, ind3) = acc_weak_av;
%     sparsity_K_lambda_mu(ind1, ind2, ind3) = mean(sum(Z ~= 0))/K(ind1);
%     tr_sparsity_K_lambda_mu(ind1, ind2, ind3) = mean(sum(Dict_mix.Z ~= 0))/K(ind1);

%     [SW,SB]=calcfisher(Dict_mix.Z,Database.tr_label,opts);
%     fWZ=trace(W'*SW*W)-trace(W'*SB*W)+norm(W'*Dict_mix.Z,'fro')^2;              
% 
%                 opts.lambda1*sum(abs(Dict_mix.Z(:)))
%                 opts.mu*fWZ
end 
end
end
end
end
% end
% save('t1_results_','result_K_lambda_mu','result_K_lambda_muWEEK',...
%     'sparsity_K_lambda_mu','tr_sparsity_K_lambda_mu')
toc

figure