
close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././FDDLPOW'))

addpath(genpath('D:\Stored_Data\data'))
SNR_INF = 2000;
cvortest = 1;  % 1 means crossvalidation data

%% testing part
% %% 
% mixture_n = 1; % mixture_n classes mixture, = 1,2,3
% pctrl.db = 0; % dynamic ratio is 0 3, 6, 10, 20 db
% pctrl.if2weak = 0; % if 2 weak components in mixture of 3
% if pctrl.db == 0     pctrl.equal = 1; else    pctrl.equal = 0; end
% Database = load_data_new(mixture_n, SNR_INF, pctrl, 1000);   
% 
% r_l1 = zeros(3,5,1); % 3 algs; L=1,2,3; 5 dynamic ratio; 5 folds
% SNR = [2000, 20, 10, 5, 0, -10, -20];
% for ind_snr = 1:length(SNR)
% %     db = Database;
% %     db.test_data = awgn(db.test_data, SNR(ind_snr), 'measured');
% db = load_data_new(mixture_n, SNR(ind_snr), pctrl, 1000);   
% for f = 1000
% for alg = 1:3
%     
%     if alg == 1 load(['FDDLOW_k100_lmbd0.0001_mu0.001_Q16.mat']);disp(opts.Dictnm); end
%     if alg == 2 
%         load('FDDLOW_mix_k100_lmbd0.0001_mu0.001_Q16_nu0.01_beta-1.mat')
%         Dict = Dict_mix;
%         disp(opts.mixnm); 
%     end
%     if alg == 3 
%         load('FDDLOW_mix_k100_lmbd0.0001_mu0.001_Q16_nu0.01_beta0.0001.mat')
%         Dict = Dict_mix;
%         disp(opts.mixnm); 
%     end
%          
%     % run prep_ZF 
%     Z = sparsecoding(Dict, db, opts, mixture_n, 0);
%     Z = aoos(Z,db.featln, size(Z, 2));
%     
%     Xtestorcv = Dict.W'*Z;
%     Xtr = Dict.W'*Dict.Z;%aoos(Dict.Z,Database.featln, size(Dict.Z, 2));
%     % KNN classifier
%     if alg == 1
%         acc = myknn(Xtr, Xtestorcv, db, 0); % k = 5 ;
%     else
%         run calc_M
%         H = W'*M;
%         r_zeroforcing = pinv(H)*W'*Z;
%         [~, labels_pre] = sort(r_zeroforcing, 1, 'descend');
%         cvtlabel = aoos(db.test_label,db.featln,size(db.test_label, 2));
%         a = labels_pre(1,:) - cvtlabel(1,:);
%         acc = length(a(a ==0))/length(a);
%     end
%     r_l1(alg, ind_snr, f-999) = acc;   
% end
% end
% end
% % save('test_L=1.mat','r_knn');
% plot(mean(r_l1,3)', '-x')



%% 
mixture_n = 2; 

r_zf = zeros(3,5,1); %3 algs; L=2; 5 dynamic ratio; 5 folds
r_mf = r_zf; r_zf_weak = r_zf; r_mf_weak = r_mf;
r_lr = r_zf; r_lr_weak = r_zf; r_nn_weak = r_zf; r_nn = r_zf;
dynamic_ratio = [0, 3, 6, 10, 20];    
for indd = 1:5
    pctrl.db = dynamic_ratio(indd); % dynamic ratio is 0 3, 6, 10, 20 db
    pctrl.if2weak = 0; % if 2 weak components in mixture of 3
    if pctrl.db == 0     pctrl.equal = 1; else    pctrl.equal = 0; end
    Database = load_data_new(mixture_n, SNR_INF, pctrl, 1000);   
for f = 1000
for alg = 1:3

    if alg == 1 
        load(['FDDLOW_k100_lmbd0.0001_mu0.001_Q16.mat']);
%         load('B_X_Y_dict1.mat')
%         load('NN_dict1.mat')
        disp(opts.Dictnm); 
    end
    if alg == 2 
        load('FDDLOW_mix_k100_lmbd0.0001_mu0.001_Q16_nu0.01_beta-1.mat')
%         load('B_X_Y_dict2.mat')
%         load('NN_dict2.mat')
        Dict = Dict_mix;
        disp(opts.mixnm);
    end
    if alg == 3 
        load('FDDLOW_mix_k100_lmbd0.0001_mu0.001_Q16_nu0.01_beta0.0001.mat')
%         load('B_X_Y_dict3.mat')
%         load('NN_dict3.mat')
        Dict = Dict_mix;
        disp(opts.mixnm); 
    end
    
    run calc_M
    Z = sparsecoding(Dict, Database, opts, mixture_n, 0);
    Z = aoos(Z,Database.featln, size(Z, 2));
    % zero forcing
    H = W'*M;
    r_zeroforcing = pinv(H)*W'*Z;
    [~, labels_pre] = sort(r_zeroforcing, 1, 'descend');

    % matched filter
    r_matched = H'*W'*Z;
    [~, labels_pre_mf] = sort(r_matched, 1, 'descend');

%     %neural networks
%     outputs = net(W'*Z);
%     [~, labels_pre_nn] = sort(outputs, 1, 'descend');
    
    % calculate accuracy
    [~, r_zf_weak(alg, indd, f-999), r_zf(alg, indd, f-999)] = calc_labels(labels_pre, opts);
    [~, r_mf_weak(alg, indd, f-999), r_mf(alg, indd, f-999)] = calc_labels(labels_pre_mf, opts);
%     [~, r_nn_weak(alg, indd, f-999), r_nn(alg, indd, f-999)] = calc_labels(labels_pre_nn, opts);
%     [~, r_lr_weak(alg, indd, f-999), r_lr(alg, indd, f-999)] = lr_test(Dict, Database, Z, B, pctrl);

end
end   
end
% save('test_L=2.mat','r_zf', 'r_zf_weak', 'r_mf', 'r_mf_weak', 'r_nn','r_nn_weak','r_lr','r_lr_weak');
figure;

%% 
mixture_n = 3;

r_zf = zeros(3,5,1); %3 algs; L=2; 5 dynamic ratio; 5 folds
r_mf = r_zf; r_zf_weak = r_zf; r_mf_weak = r_mf;
r_lr = r_zf; r_lr_weak = r_zf; r_nn_weak = r_zf; r_nn = r_zf;
dynamic_ratio = [0, 3, 6, 10, 20];    
for indd = 1:5
    pctrl.db = dynamic_ratio(indd); % dynamic ratio is 0 3, 6, 10, 20 db
    pctrl.if2weak = 0; % if 2 weak components in mixture of 3
    if pctrl.db == 0     pctrl.equal = 1; else    pctrl.equal = 0; end
    Database = load_data_new(mixture_n, SNR_INF, pctrl, 1000);   
for f = 1000
for alg = 1:3

    if alg == 1 
        load(['FDDLOW_k100_lmbd0.0001_mu0.001_Q16.mat']);
%         load('B_X_Y_dict1.mat')
%         load('NN_dict1.mat')
        disp(opts.Dictnm); 
    end
    if alg == 2 
        load('FDDLOW_mix_k100_lmbd0.0001_mu0.001_Q16_nu0.01_beta-1.mat')
%         load('B_X_Y_dict2.mat')
%         load('NN_dict2.mat')
        Dict = Dict_mix;
        disp(opts.mixnm);
    end
    if alg == 3 
        load('FDDLOW_mix_k100_lmbd0.0001_mu0.001_Q16_nu0.01_beta0.0001.mat')
%         load('B_X_Y_dict3.mat')
%         load('NN_dict3.mat')
        Dict = Dict_mix;
        disp(opts.mixnm); 
    end
    
    run calc_M
    Z = sparsecoding(Dict, Database, opts, mixture_n, 0);
    Z = aoos(Z,Database.featln, size(Z, 2));
    % zero forcing
    H = W'*M;
    r_zeroforcing = pinv(H)*W'*Z;
    [~, labels_pre] = sort(r_zeroforcing, 1, 'descend');

    % matched filter
    r_matched = H'*W'*Z;
    [~, labels_pre_mf] = sort(r_matched, 1, 'descend');

%     %neural networks
%     outputs = net(W'*Z);
%     [~, labels_pre_nn] = sort(outputs, 1, 'descend');
    
    % calculate accuracy
    [~, r_zf_weak(alg, indd, f-999), r_zf(alg, indd, f-999)] = calc_labels(labels_pre, opts);
    [~, r_mf_weak(alg, indd, f-999), r_mf(alg, indd, f-999)] = calc_labels(labels_pre_mf, opts);
%     [~, r_nn_weak(alg, indd, f-999), r_nn(alg, indd, f-999)] = calc_labels(labels_pre_nn, opts);
%     [~, r_lr_weak(alg, indd, f-999), r_lr(alg, indd, f-999)] = lr_test(Dict, Database, Z, B, pctrl);

end
end   
end
% save('test_L=3.mat','r_zf', 'r_zf_weak', 'r_mf', 'r_mf_weak', 'r_nn','r_nn_weak','r_lr','r_lr_weak');
figure;


