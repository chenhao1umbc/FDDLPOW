% {
close all
clear
clc;
tic
addpath(genpath('/extra/chenhao1/FDDLOW/fddlow'))
addpath(genpath('/extra/chenhao1/FDDLOW/data'))
addpath(genpath('/extra/chenhao1/DICTOL-master'))


% do traing or do crossvalidation
cvortest = 0;  % 1 means cv, 0 means test
SNR_INF = 2000;
K = 25;
T0=  15;
C=6;
param.L = T0; % not more than 10 non-zeros coefficients
param.eps=0.0; % squared norm of the residual should be less than 0.1
load('NN_ksvd.mat')
load('B_X_Y_ksvd.mat')
zf.acc_week = zeros(5,5);
zf.acc = zeros(5,5);
mf.acc_week = zeros(5,5);
mf.acc = zeros(5,5);
nn.acc_week = zeros(5,5);
nn.acc = zeros(5,5);
lr.acc_week = zeros(5,5);
lr.acc = zeros(5,5);

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
    
    wr_zf = pinv(H)*Z;
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


%%%%%%%%%%%plot results%%%%%%%%%%%%%%%%%%%

zfw = mean(zf.acc_week, 2);
zfa = mean(zf.acc, 2);
mfa = mean(mf.acc, 2);
mfw = mean(mf.acc_week, 2);
nna = mean(nn.acc, 2);
nnw = mean(nn.acc_week, 2);
lra = mean(lr.acc, 2);
lrw = mean(lr.acc_week, 2);

figure;
plot([0,3,6,10,20], mfa, '-x')
hold on
plot([0,3,6,10,20], mfw, '--o')
ylim([0.2,1])
xlabel('power ratio 0, 3,6, 10, 20db')
ylabel('accuracy')
legend('all', 'weak')
title(['L=', num2str(mixture_n), 'ksvd + mf'])

figure;
plot([0,3,6,10,20], zfa, '-x')
hold on
plot([0,3,6,10,20], zfw, '--o')
ylim([0.2,1])
xlabel('power ratio 0, 3,6, 10, 20db')
ylabel('accuracy')
legend('all', 'weak')
title(['L=', num2str(mixture_n), 'ksvd + zf'])

figure;
plot([0,3,6,10,20], nna, '-x')
hold on
plot([0,3,6,10,20], nnw, '--o')
ylim([0.2,1])
xlabel('power ratio 0, 3,6, 10, 20db')
ylabel('accuracy')
legend('all', 'weak')
title(['L=', num2str(mixture_n), 'ksvd + nn'])

figure;
plot([0,3,6,10,20], lra, '-x')
hold on
plot([0,3,6,10,20], lrw, '--o')
ylim([0.2,1])
xlabel('power ratio 0, 3,6, 10, 20db')
ylabel('accuracy')
legend('all', 'weak')
title(['L=', num2str(mixture_n), 'ksvd + lr'])

%}

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script is used to test/cv LRSDL+SVM for mixture=1
close all
clear
clc;
tic

addpath(genpath('/extra/chenhao1/FDDLOW/fddlow'))
addpath(genpath('/extra/chenhao1/FDDLOW/data'))
addpath(genpath('/extra/chenhao1/DICTOL-master'))

%% load data
mixture_n = 1; % mixture_n classes mixture, = 1,2,3
SNR_INF = 2000;
pctrl.db = 20; % dynamic ratio is 0 3, 6, 10, 20 db
if mixture_n < 3  pctrl.if2weak = 0; end
if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end

for f = 1000:1004
[Database]=load_data_new(mixture_n, SNR_INF, pctrl, f);

%% settings
K = 25;
T0=  18;

params = ['k_',num2str(K), 'T0_',num2str(T0), 'f_', num2str(f)];
load([params,'ksvd_train.mat'])


%% test/cv
X = Database.cv_data;
param.L = T0; % not more than 10 non-zeros coefficients
param.eps=0.0; % squared norm of the residual should be less than 0.1
Z = mexOMP(X,D,param);
Z = aoos(Z,Database.featln, size(Z, 2)); 
acc_knn(f-999) = myknn(W, Z, Database, 0) % k = 5 ;

end % end of f=1000:1004
sum(acc_knn)/5
% end  % end of i 
%}