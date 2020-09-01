% draft file for testing codes

close all
clear
clc;
tic

addpath(genpath('../../fddlow'))
addpath(genpath('../../data'))
addpath(genpath('../../FDDLPOW'))

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
title(['L=', num2str(mixture_n), 'alg123 + svm'])

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
title(['L=', num2str(mixture_n), 'alg123 + svm'])

%% _____________________KSVD svm_____________________

close all
clear
clc;
tic
addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././DICTOL-master'))

% do traing or do crossvalidation
cvortest = 0;  % 1 means cv, 0 means test
SNR_INF = 2000;
K = 25;
T0=  15;
C=6;
param.L = T0; % not more than 10 non-zeros coefficients
param.eps=0.0; % squared norm of the residual should be less than 0.1
load('Mdl_X_Y_ksvd.mat')
svm.acc_week = zeros(5,5);
svm.acc = zeros(5,5);

mixture_n = 2; % mixture_n classes mixture, = 1,2,3
dbpool = [0 3 6 10 20];
for ind = 1:5
pctrl.db = dbpool(ind); % dynamic ratio is 0 3, 6, 10, 20 db
if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end
if mixture_n < 3  pctrl.if2weak = 0; end
if mixture_n == 3 pctrl.if2weak = 0; end

for f = 1000:1004     
    [Database]=load_data_new(mixture_n, SNR_INF, pctrl, f);
    if cvortest== 1
        Y = Database.cv_mixdata;
    else
        Y = Database.test_mixdata;
    end
    params = ['k_',num2str(K), 'T0_',num2str(T0), 'f_', num2str(f)];
    load([params,'ksvd_train.mat'])
    X = W;
    
    % run prep_ZF 
    Z = mexOMP(Y,D,param);
    Z = aoos(Z,Database.featln, size(Z, 2));      
    
    N = size(Database.tr_label,2);
    Nc = N / C;
    opts.C = C; % 6 classes
    featln = Database.featln;
    opts.n = Database.N_c;                
    % H0 = kron(eye(C),ones(Nc, 1)/Nc); % M = Z*H3
    H0 = kron(eye(6),ones(2400, 1)/2400); 
    H = X*H0;
    if cvortest == 1
        N_t = size(Database.cv_mixlabel, 2); %%%%% cv or test **************
    else 
        N_t = size(Database.test_mixlabel, 2); %%%%% cv or test **************
    end
    opts.ln_test = N_t/featln;
    opts.equal = pctrl.equal;
    opts.Ncombs = max(Database.cv_mixlabel);
    
    [~, pre_prob, ~] = predict(Mdl, Z');
    [~,labels_pre] = sort(pre_prob, 2, 'descend');
    labels_pre = labels_pre';
    
    if pctrl.if2weak == 0
        [~, svm.acc_week(ind, f-999), svm.acc(ind, f-999)] = calc_labels(labels_pre, opts);
    else
        [~, svm.acc_week(ind, f-999), svm.acc(ind, f-999)] = calc_labels2w(labels_pre, opts);
    end

end
svm.acc
svm.acc_week
end
toc


%%%%%%%%%%%plot results%%%%%%%%%%%%%%%%%%%
zfa = mean(svm.acc, 2);
zfw = mean(svm.acc_week, 2);

figure;
plot([0,3,6,10,20], zfa, '-x')
hold on
plot([0,3,6,10,20], zfw, '--o')
ylim([0.2,1])
xlabel('power ratio 0, 3,6, 10, 20db')
ylabel('accuracy')
legend('all', 'weak')
title(['L=', num2str(mixture_n), 'ksvd + svm'])
savefig(['L=', num2str(mixture_n), 'ksvd + svm.fig'])

%%
% do traing or do crossvalidation
cvortest = 0;  % 1 means cv, 0 means test
SNR_INF = 2000;
K = 25;
T0=  15;
C=6;
param.L = T0; % not more than 10 non-zeros coefficients
param.eps=0.0; % squared norm of the residual should be less than 0.1
load('Mdl_X_Y_ksvd.mat')
svm.acc_week = zeros(5,5);
svm.acc = zeros(5,5);

mixture_n = 3; % mixture_n classes mixture, = 1,2,3
dbpool = [0 3 6 10 20];
for ind = 1:5
pctrl.db = dbpool(ind); % dynamic ratio is 0 3, 6, 10, 20 db
if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end
if mixture_n < 3  pctrl.if2weak = 0; end
if mixture_n == 3 pctrl.if2weak = 0; end

for f = 1000:1004     
    [Database]=load_data_new(mixture_n, SNR_INF, pctrl, f);
    if cvortest== 1
        Y = Database.cv_mixdata;
    else
        Y = Database.test_mixdata;
    end
    params = ['k_',num2str(K), 'T0_',num2str(T0), 'f_', num2str(f)];
    load([params,'ksvd_train.mat'])
    X = W;
    
    % run prep_ZF 
    Z = mexOMP(Y,D,param);
    Z = aoos(Z,Database.featln, size(Z, 2));      
    
    N = size(Database.tr_label,2);
    Nc = N / C;
    opts.C = C; % 6 classes
    featln = Database.featln;
    opts.n = Database.N_c;                
    % H0 = kron(eye(C),ones(Nc, 1)/Nc); % M = Z*H3
    H0 = kron(eye(6),ones(2400, 1)/2400); 
    H = X*H0;
    if cvortest == 1
        N_t = size(Database.cv_mixlabel, 2); %%%%% cv or test **************
    else 
        N_t = size(Database.test_mixlabel, 2); %%%%% cv or test **************
    end
    opts.ln_test = N_t/featln;
    opts.equal = pctrl.equal;
    opts.Ncombs = max(Database.cv_mixlabel);
    
    [~, pre_prob, ~] = predict(Mdl, Z');
    [~,labels_pre] = sort(pre_prob, 2, 'descend');
    labels_pre = labels_pre';
    
    if pctrl.if2weak == 0
        [~, svm.acc_week(ind, f-999), svm.acc(ind, f-999)] = calc_labels(labels_pre, opts);
    else
        [~, svm.acc_week(ind, f-999), svm.acc(ind, f-999)] = calc_labels2w(labels_pre, opts);
    end

end
svm.acc
svm.acc_week
end
toc


%%%%%%%%%%%plot results%%%%%%%%%%%%%%%%%%%
zfa = mean(svm.acc, 2);
zfw = mean(svm.acc_week, 2);

figure;
plot([0,3,6,10,20], zfa, '-x')
hold on
plot([0,3,6,10,20], zfw, '--o')
ylim([0.2,1])
xlabel('power ratio 0, 3,6, 10, 20db')
ylabel('accuracy')
legend('all', 'weak')
title(['L=', num2str(mixture_n), 'ksvd + svm'])

savefig(['L=', num2str(mixture_n), 'ksvd + svm.fig'])


%% _________________________lrsdl svm___________________
% {
close all
clear
clc;
tic
addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././DICTOL-master'))

% do traing or do crossvalidation
cvortest = 0;  % 1 means cv, 0 means test
SNR_INF = 2000;
lbmd_bar = 0.001;
k0 = 3;
k = 4;
lambda1 = 1e-4;
lambda2 = 5e-3;
lambda3 = 5e-2;
C = 6;
load('Mdl_X_Y_lrsdl.mat')
svm.acc_week = zeros(5,5);
svm.acc = zeros(5,5);


%% testing/cv part
% load data
mixture_n = 2; % mixture_n classes mixture, = 1,2,3

dbpool = [0 3 6 10 20];
for ind = 1:5
pctrl.db = dbpool(ind); % dynamic ratio is 0 3, 6, 10, 20 db
if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end
if mixture_n < 3  pctrl.if2weak = 0; end
if mixture_n == 3 pctrl.if2weak = 0; end

for f = 1000:1004     
    [Database]=load_data_new(mixture_n, SNR_INF, pctrl, f);
    if cvortest== 1
        Y = Database.cv_mixdata;
    else
        Y = Database.test_mixdata;
    end
    param = ['k_',num2str(k),'k0_',num2str(k0), 'l1_',num2str(lambda1), ...
        'l2_',num2str(lambda2), 'l3_',num2str(lambda3), 'f_', num2str(f)];
    filename = [param,'lrscdl_train.mat'];
    load(filename)    
    
    % run prep_ZF 
    [Z, Z0] = local_sparse_coding(Y, D, D0, CoefM0, lambda1, lambda2);
    Z = aoos(Z,Database.featln, size(Z, 2));      
    
    N = size(Database.tr_label,2);
    Nc = N / C;
    opts.C = C; % 6 classes
    featln = Database.featln;
    opts.n = Database.N_c;                
    % H0 = kron(eye(C),ones(Nc, 1)/Nc); % M = Z*H3
    H0 = kron(eye(6),ones(2400, 1)/2400); 
    H = X*H0;
    if cvortest == 1
        N_t = size(Database.cv_mixlabel, 2); %%%%% cv or test **************
    else 
        N_t = size(Database.test_mixlabel, 2); %%%%% cv or test **************
    end
    opts.ln_test = N_t/featln;
    opts.equal = pctrl.equal;
    opts.Ncombs = max(Database.cv_mixlabel);
       
    [~, pre_prob, ~] = predict(Mdl, Z');
    [~,labels_pre] = sort(pre_prob, 2, 'descend');
    labels_pre = labels_pre';
    
    if pctrl.if2weak == 0
        [~, svm.acc_week(ind, f-999), svm.acc(ind, f-999)] = calc_labels(labels_pre, opts);
    else
        [~, svm.acc_week(ind, f-999), svm.acc(ind, f-999)] = calc_labels2w(labels_pre, opts);
    end
end
svm.acc
svm.acc_week

end
toc


%%%%%%%%%%%plot results%%%%%%%%%%%%%%%%%%%

zfw = mean(svm.acc_week, 2);
zfa = mean(svm.acc, 2);


figure;
plot([0,3,6,10,20], zfa, '-x')
hold on
plot([0,3,6,10,20], zfw, '--o')
ylim([0.2,1])
xlabel('power ratio 0, 3,6, 10, 20db')
ylabel('accuracy')
legend('all', 'weak')
title(['L=', num2str(mixture_n), 'lrsdl + svm'])
savefig(['L=', num2str(mixture_n), 'lrsdl + svm.fig'])

%}

%% 
svm.acc_week = zeros(5,5);
svm.acc = zeros(5,5);


%% testing/cv part
% load data
mixture_n = 3; % mixture_n classes mixture, = 1,2,3

dbpool = [0 3 6 10 20];
for ind = 1:5
pctrl.db = dbpool(ind); % dynamic ratio is 0 3, 6, 10, 20 db
if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end
if mixture_n < 3  pctrl.if2weak = 0; end
if mixture_n == 3 pctrl.if2weak = 0; end

for f = 1000:1004     
    [Database]=load_data_new(mixture_n, SNR_INF, pctrl, f);
    if cvortest== 1
        Y = Database.cv_mixdata;
    else
        Y = Database.test_mixdata;
    end
    param = ['k_',num2str(k),'k0_',num2str(k0), 'l1_',num2str(lambda1), ...
        'l2_',num2str(lambda2), 'l3_',num2str(lambda3), 'f_', num2str(f)];
    filename = [param,'lrscdl_train.mat'];
    load(filename)    
    
    % run prep_ZF 
    [Z, Z0] = local_sparse_coding(Y, D, D0, CoefM0, lambda1, lambda2);
    Z = aoos(Z,Database.featln, size(Z, 2));      
    
    N = size(Database.tr_label,2);
    Nc = N / C;
    opts.C = C; % 6 classes
    featln = Database.featln;
    opts.n = Database.N_c;                
    % H0 = kron(eye(C),ones(Nc, 1)/Nc); % M = Z*H3
    H0 = kron(eye(6),ones(2400, 1)/2400); 
    H = X*H0;
    if cvortest == 1
        N_t = size(Database.cv_mixlabel, 2); %%%%% cv or test **************
    else 
        N_t = size(Database.test_mixlabel, 2); %%%%% cv or test **************
    end
    opts.ln_test = N_t/featln;
    opts.equal = pctrl.equal;
    opts.Ncombs = max(Database.cv_mixlabel);
       
    [~, pre_prob, ~] = predict(Mdl, Z');
    [~,labels_pre] = sort(pre_prob, 2, 'descend');
    labels_pre = labels_pre';
    
    if pctrl.if2weak == 0
        [~, svm.acc_week(ind, f-999), svm.acc(ind, f-999)] = calc_labels(labels_pre, opts);
    else
        [~, svm.acc_week(ind, f-999), svm.acc(ind, f-999)] = calc_labels2w(labels_pre, opts);
    end
end
svm.acc
svm.acc_week

end
toc


%%%%%%%%%%%plot results%%%%%%%%%%%%%%%%%%%

zfw = mean(svm.acc_week, 2);
zfa = mean(svm.acc, 2);


figure;
plot([0,3,6,10,20], zfa, '-x')
hold on
plot([0,3,6,10,20], zfw, '--o')
ylim([0.2,1])
xlabel('power ratio 0, 3,6, 10, 20db')
ylabel('accuracy')
legend('all', 'weak')
title(['L=', num2str(mixture_n), 'lrsdl + svm'])
savefig(['L=', num2str(mixture_n), 'lrsdl + svm.fig'])