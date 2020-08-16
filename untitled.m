% draft file for testing codes

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

mixture_n = 2; 

r_svm = zeros(3,5,5); %3 algs; L=2; 5 dynamic ratio; 5 folds
r_svm_weak = r_svm; 

dynamic_ratio = [0, 3, 6, 10, 20];    
for indd = 1:5
    pctrl.db = dynamic_ratio(indd); % dynamic ratio is 0 3, 6, 10, 20 db
    pctrl.if2weak = 0; % if 2 weak components in mixture of 3
    if pctrl.db == 0     pctrl.equal = 1; else    pctrl.equal = 0; end
    Database = load_data_new(mixture_n, SNR_INF, pctrl, 1000);   
for f = 1000:1004
for alg = 1:3

    if alg == 1 
        load(['dict1_k25_lmbd0.01_mu0.1_Q10_rng',num2str(f),'.mat']);
        load('Mdl_X_Y_alg1.mat')
        opts.lambda1 = 0.01/8;
        disp(opts.Dictnm); 
    end
    if alg == 2 
        load(['dict2_k25_lmbd0.1_mu0.001_Q20_nu10_rng',num2str(f+5),'.mat']);
        load('Mdl_X_Y_alg2.mat')
        opts.lambda1 = 0.025;
        disp(opts.Dict2nm); 
    end
    if alg == 3 
        load(['dict3_k25_lmbd0.1_mu0.001_Q20_nu10_beta1_rng',num2str(f+5),'.mat']);
        load('Mdl_X_Y_alg3.mat')
        opts.lambda1 = 0.025;
        disp(opts.Dict3nm); 
    end
    
    run calc_M
    Z = sparsecoding(Dict, Database, opts, mixture_n, 0);
    Z = aoos(Z,Database.featln, size(Z, 2));

    wz = W'*Z; %(W'*aoos(Z, featln, N));
    [~, pre_prob, ~] = predict(Mdl, wz');
    [~,labels_pre] = sort(pre_prob, 2, 'descend');
    labels_pre = labels_pre';

    if pctrl.if2weak == 0
        [acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts);
    else
        [acc_weak, acc_weak_av, acc_all] = calc_labels2w(labels_pre, opts);
    end
    r_svm(alg, indd, f-999)= acc_all;
    r_svm_weak(alg, indd, f-999)= acc_weak_av;

end
end   
end
save(['test_L=',num2str(mixture_n),'.mat'],'r_svm', 'r_svm_weak');


r_svm_weak = mean(r_svm_weak, 3);
r_svm = mean(r_svm, 3);

figure;
for alg = 1:3
plot([0,3,6,10,20], r_svm(alg,:), '-x')
hold on
plot([0,3,6,10,20], r_svm_weak(alg,:), '--o')
ylim([0.2,1])
xlabel('power ratio 0, 3,6, 10, 20db')
ylabel('accuracy')
end
legend('all', 'weak')
title(['L=', num2str(mixture_n), 'alg124 + svm'])
%%
mixture_n = 3; 

r_svm = zeros(3,5,5); %3 algs; L=2; 5 dynamic ratio; 5 folds
r_svm_weak = r_svm; 

dynamic_ratio = [0, 3, 6, 10, 20];    
for indd = 1:5
    pctrl.db = dynamic_ratio(indd); % dynamic ratio is 0 3, 6, 10, 20 db
    pctrl.if2weak = 0; % if 2 weak components in mixture of 3
    if pctrl.db == 0     pctrl.equal = 1; else    pctrl.equal = 0; end
    Database = load_data_new(mixture_n, SNR_INF, pctrl, 1000);   
for f = 1000:1004
for alg = 1:3

    if alg == 1 
        load(['dict1_k25_lmbd0.01_mu0.1_Q10_rng',num2str(f),'.mat']);
        load('Mdl_X_Y_alg1.mat')
        opts.lambda1 = 0.01/8;
        disp(opts.Dictnm); 
    end
    if alg == 2 
        load(['dict2_k25_lmbd0.1_mu0.001_Q20_nu10_rng',num2str(f),'.mat']);
        load('Mdl_X_Y_alg2.mat')
        opts.lambda1 = 0.025;
        disp(opts.Dict2nm); 
    end
    if alg == 3 
        load(['dict3_k25_lmbd0.1_mu0.001_Q20_nu10_beta1_rng',num2str(f+5),'.mat']);
        load('Mdl_X_Y_alg3.mat')
        opts.lambda1 = 0.025;
        disp(opts.Dict3nm); 
    end
    
    run calc_M
    Z = sparsecoding(Dict, Database, opts, mixture_n, 0);
    Z = aoos(Z,Database.featln, size(Z, 2));

    wz = W'*Z; %(W'*aoos(Z, featln, N));
    [~, pre_prob, ~] = predict(Mdl, wz');
    [~,labels_pre] = sort(pre_prob, 2, 'descend');
    labels_pre = labels_pre';

    if pctrl.if2weak == 0
        [acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts);
    else
        [acc_weak, acc_weak_av, acc_all] = calc_labels2w(labels_pre, opts);
    end
    r_svm(alg, indd, f-999)= acc_all;
    r_svm_weak(alg, indd, f-999)= acc_weak_av;

end
end   
end
save(['test_L=',num2str(mixture_n),'.mat'],'r_svm', 'r_svm_weak');

r_svm_weak = mean(r_svm_weak, 3);
r_svm = mean(r_svm, 3);


figure;
for alg = 1:3
plot([0,3,6,10,20], r_svm(alg,:), '-x')
hold on
plot([0,3,6,10,20], r_svm_weak(alg,:), '--o')
ylim([0.2,1])
xlabel('power ratio 0, 3,6, 10, 20db')
ylabel('accuracy')
end
legend('all', 'weak')
title(['L=', num2str(3), 'alg124 + svm'])