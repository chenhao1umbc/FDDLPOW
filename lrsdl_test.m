% draft file for testing codes

close all
clear
clc;
tic
addpath(genpath('/extra/chenhao1/FDDLOW/fddlow'))
addpath(genpath('/extra/chenhao1/FDDLOW/data'))
addpath(genpath('/extra/chenhao1/DICTOL-master'))

% do traing or do crossvalidation
cvortest = 1;  % 1 means cv, 0 means test
SNR_INF = 2000;
lbmd_bar = 0.001;
k0 = 3;
k = 15;
lambda1 = 1e-4;
lambda2 = 1e-3;
lambda3 = 1e-2;
C = 6;
load('NN_lrsdl.mat')
load('B_X_Y_lrsdl.mat')
zf.acc_week = zeros(5,5);
zf.acc = zeros(5,5);
mf.acc_week = zeros(5,5);
mf.acc = zeros(5,5);
nn.acc_week = zeros(5,5);
nn.acc = zeros(5,5);
lr.acc_week = zeros(5,5);
lr.acc = zeros(5,5);

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
    [Z, Z0] = local_sparse_coding(Y, D, D0, CoefM0, lbmd_bar, lambda2);
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
    
    r_zf = pinv(H)*Z;
    [~, labels_pre_zf] = sort(r_zf, 1, 'descend');
    % matched filter
    r_matched = H'*Z;
    [~, labels_pre_mf] = sort(r_matched, 1, 'descend');
    %neural networks
    outputs = net(Z);
    [~, labels_pre_nn] = sort(outputs, 1, 'descend');
    
    [~, zf.acc_week(ind, f-999), zf.acc(ind, f-999)] = calc_labels(labels_pre_zf, opts);
    [~, mf.acc_week(ind, f-999), mf.acc(ind, f-999)] = calc_labels(labels_pre_mf, opts);
    [~, nn.acc_week(ind, f-999), nn.acc(ind, f-999)] = calc_labels(labels_pre_nn, opts);
    
    pre_prob = mnrval(B, Z');
    [~,labels_pre] = sort(pre_prob, 2, 'descend');
    labels_pre = labels_pre';
    if pctrl.if2weak == 0
        [~, lr.acc_week(ind, f-999), lr.acc(ind, f-999)] = calc_labels(labels_pre, opts);
    else
        [~, lr.acc_week(ind, f-999), lr.acc(ind, f-999)] = calc_labels2w(labels_pre, opts);
    end

end
zf.acc
zf.acc_week
mf.acc
mf.acc_week
nn.acc
nn.acc_week
lr.acc
lr.acc_week
end
toc
