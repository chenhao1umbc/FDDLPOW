
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

% load settings
K = [25, 50, 100, 150];
lbmd = [0.1, 0.01, 0.001, 00001];
mu = [0.1, 0.01, 0.01, 0.001, 0.0001];
% SNR = [2000, 20, 0, -5, -10, -20];
Q = [10 20 25 30 40 50 75 100 150];

SNR = 20;
r = zeros(length(mu),length(lbmd), length(Q), length(K), 5); % Q, lambda, folds 
r_zf = r; r_mf = r;
%% CV/testing part
for f = 1000:1004
[Database]=load_data_new(mixture_n, SNR, pctrl, f);

if do_cv ==1      
% result_K_lambda_mu = zeros(length(K),length(lbmd),length(mu));
% sparsity_K_lambda_mu = zeros(length(K),length(lbmd),length(mu));
% tr_sparsity_K_lambda_mu = zeros(length(K),length(lbmd),length(mu));
% result_K_lambda_muWEEK = zeros(length(K),length(lbmd),length(mu));

for indk = 1: length(K)
for indq = 1: length(Q) 
for indl = 1: length(lbmd)
for indm = 1: length(mu)

    [opts]=loadoptions(K(indk),lbmd(indl),mu(indm),Q(indq),1000,1, SNR, f);
    disp(opts.Dictnm)
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
        acc_knn = myknn(Xtr, Xtestorcv, Database, cvortest) % k = 5 ;
        r(indm, indl, indq, indk, f-999) = acc_knn;
    end

    run calc_M
    % zero forcing
    H = W'*M;
    r_zf = pinv(H)*W'*Z;
    [~, labels_pre] = sort(r_zf, 1, 'descend');
    
    % matched filter
    r_matched = H'*W'*Z;
    [~, labels_pre_mf] = sort(r_zf, 1, 'descend');
    
    % calculate accuracy
    if Database.N_c == 1
        acc_ZF = myZFMF(labels_pre, Database, cvortest)
        acc_MF = myZFMF(labels_pre_mf, Database, cvortest)
        r_zf(indm, indl, indq, indk, f-999) = acc_ZF;
        r_mf(indm, indl, indq, indk, f-999) = acc_MF;
    else
        [acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts);
        [t, tt, ttt] = calc_labels(labels_pre_mf, opts);
    end
    
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
end
% save('t1_results_','result_K_lambda_mu','result_K_lambda_muWEEK',...
%     'sparsity_K_lambda_mu','tr_sparsity_K_lambda_mu')
toc

figure


% % table 3
% 
% close all
% clear
% clc;
% tic
% 
% addpath(genpath('.././fddlow'))
% addpath(genpath('.././data'))
% addpath(genpath('.././FDDLPOW'))
% 
% cv = 0; % validation or testing
% if2weak = 1; % two weak components in the mixture
% Alg_n = 3;
% for mixn = 3%[2, 3]
% for id = [3, 6, 10, 20]
% % load data
% L = mixn; % mixture_n classes mixture, = 1,2,3
% SNR = 2000;
% pctrl.db = id; % dynamic ratio is 0 3, 6, 10, 20 db
% pctrl.if2weak = if2weak; % only works for mixture_n == 3
% if pctrl.db == 0
%     pctrl.equal = 1;
% else
%     pctrl.equal = 0;
% end
% % the equal power mixture, 400 samples per combination
% [Database]=load_data_new(L, SNR, pctrl);
% 
% %% load settings
% K = 100;
% lbmd = 1e-4;
% mu=1e-3;
% Q=16;% this is wq without negative
% SNR = 2000;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % from table one we know that there are combinationes accuracy is above 0.99
% % one is K = 100, lambda = 1e-4, mu = 1e-3, nu = 0.01
% % another is K = 100, lambda = 1e-3, mu = 0.1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nu= 0.01 ;
% if Alg_n == 2
%     beta = -1;
% else
%     beta = 0.0001; % beta = -1, for alg2
% end
% 
% [opts]=loadoptions(K,lbmd,mu,Q,nu,beta, SNR);
% if SNR == 2000   
% if Alg_n ==1
%     load(opts.Dictnm)
% else
%     load(opts.mixnm) 
% end
% else
% nm = ['SNR', num2str(SNR), opts.mixnm];
% load(nm)
% end
%                            
% %% run prep_ZF 
% if exist('Dict')==1
%     Dict_mix = Dict;
% end
% if cv == 1
%     Z = sparsecoding_mix_cv(Dict_mix, Database, opts); %%%%% cv or test **************
% else
%     Z = sparsecoding_mix_test(Dict_mix, Database, opts);
% end
% W = Dict_mix.W;
% C = max(Database.tr_label);
% N = size(Database.tr_label,2);
% Nc = N / C;
% opts.C = C; % 6 classes
% featln = Database.featln;
% opts.n = Database.N_c;                
% H3 = kron(eye(C),ones(Nc, 1)/Nc); % M = Z*H3
% M = Dict_mix.Z*H3;
% % zero forcing
% H = W'*M;
% result = pinv(H)*W'*aoos(Z,featln,size(Z, 2));
% [~, labels_pre] = sort(result, 1, 'descend');
% 
% opts.Ncombs = max(Database.cv_mixlabel);
% if cv == 1
%     N_t = size(Database.cv_mixlabel, 2); %%%%% cv or test **************
% else 
%     N_t = size(Database.test_mixlabel, 2); %%%%% cv or test **************
% end
% opts.ln_test = N_t/featln;
% opts.equal = pctrl.equal;
% 
% if if2weak ==1
% [acc_weak, acc_weak_av, acc_all] = calc_labels2w(labels_pre, opts);  
% else
% [acc_strong, acc_weak_av, acc_all] = calc_labels(labels_pre, opts);  
% end
% acc_all
% acc_weak_av
% end
% end
% toc
