
close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././FDDLPOW'))

addpath(genpath('D:\Stored_Data\data'))
SNR_INF = 2000;
cvortest = 0;  % 1 means crossvalidation data

%% testing part


%% 
mixture_n = 3; 

r_zf = zeros(3,5,5); %3 algs; L=2; 5 dynamic ratio; 5 folds
r_mf = r_zf; r_zf_weak = r_zf; r_mf_weak = r_mf;
r_lr = r_zf; r_lr_weak = r_zf; r_nn_weak = r_zf; r_nn = r_zf;
dynamic_ratio = [3, 6, 10, 20];    
for indd = 1:4
    pctrl.db = dynamic_ratio(indd); % dynamic ratio is 0 3, 6, 10, 20 db
    pctrl.if2weak = 1; % if 2 weak components in mixture of 3
    if pctrl.db == 0     pctrl.equal = 1; else    pctrl.equal = 0; end
    Database = load_data_new(mixture_n, SNR_INF, pctrl, 1000);   
for f = 1000:1004
for alg = 1:3

    if alg == 1 
        load(['dict1_k25_lmbd0.01_mu0.1_Q10_rng',num2str(f),'.mat']);
        load('B_X_Y_dict1.mat')
        load('NN_dict1.mat')
        opts.lambda1 = 0.01/8;
        disp(opts.Dictnm); 
    end
    if alg == 2 
        load(['dict2_k25_lmbd0.1_mu0.001_Q20_nu10_rng',num2str(f),'.mat']);
        load('B_X_Y_dict2.mat')
        load('NN_dict2.mat')
        opts.lambda1 = 0.05;
        disp(opts.Dict2nm); 
    end
    if alg == 3 
        load(['dict3_k25_lmbd0.1_mu0.001_Q20_nu10_beta1_rng',num2str(f),'.mat']);
        load('B_X_Y_dict3.mat')
        load('NN_dict3.mat')
        opts.lambda1 = 0.025;
        disp(opts.Dict3nm); 
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

    %neural networks
    outputs = net(W'*Z);
    [~, labels_pre_nn] = sort(outputs, 1, 'descend');
    
    % calculate accuracy
    [~, r_zf_weak(alg, indd, f-999), r_zf(alg, indd, f-999)] = calc_labels2w(labels_pre, opts);
    [~, r_mf_weak(alg, indd, f-999), r_mf(alg, indd, f-999)] = calc_labels2w(labels_pre_mf, opts);
    [~, r_nn_weak(alg, indd, f-999), r_nn(alg, indd, f-999)] = calc_labels2w(labels_pre_nn, opts);
    [~, r_lr_weak(alg, indd, f-999), r_lr(alg, indd, f-999)] = lr_test(Dict, Database, Z, B, pctrl);

end
end   
end
save('2w_test_L=3.mat','r_zf', 'r_zf_weak', 'r_mf', 'r_mf_weak', 'r_nn','r_nn_weak','r_lr','r_lr_weak');



