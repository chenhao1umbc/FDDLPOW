% table 1, 2,3 test data only no training

close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././FDDLPOW'))

%% load settings
% do traing or do crossvalidation
cvortest = [0, 1]; % [docv, dotest] cannot be [1, 1]
alg_n = 3; % algorithm number
pctrl.db = 15; % dynamic ratio is 0 3, 5, 10, 15 db for L>=2
L = 1; % mixture_n classes mixture, = 1,2 (1 means non -mixture)
if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end
K = 60;
lbmd = 0.025;
mu= 0.005;
Q = 0.9;% this is wq without negative
nu = 0.03;
beta = 0.01;
SNR = 2000;

%% load data
[Database] = load_ESC(L, SNR, pctrl);

%% cross-val part
addpath(genpath('.././tempresult'))
[opts]=loadoptions_ESC(alg_n ,K,lbmd,mu,Q*K,nu, beta);    
if exist(opts.Dictnm, 'file')
load(opts.Dictnm,'Dict','opts')
Z = sparsecoding(Dict,Database,opts,L, cvortest);
Z = aoos(Z,Database.featln, size(Z, 2));
Xtestorcv = Dict.W'*Z;
Xtr = Dict.W'*Dict.Z;

W = Dict.W;
C = max(Database.tr_label);
N = size(Database.tr_label,2);
Nc = N / C;
opts.C = C; % 10 classes
featln = Database.featln;
opts.n = Database.N_c;                
H3 = kron(eye(C),ones(Nc, 1)/Nc); % M = Z*H3
M = Dict.Z*H3;
% zero forcing
H = W'*M;
result = pinv(H)*W'*Z;
[~, labels_pre] = sort(result, 1, 'descend');
opts.Ncombs = max(Database.cv_mixlabel);
opts.ln_test = size(Database.test_mixlabel, 2)/featln;
opts.equal = pctrl.equal;
% run zero-forcing
[acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts);
% run mlknn
[acc_weak_mlsvm, acc_weak_av_mlsvm, acc_all_mlsvm] = mymlsvm(aoos(Xtr,...
    Database.featln, size(Xtr, 2)), Xtestorcv, cvortest, opts);
% % run mlknn
% [acc_weak_mlknn, acc_weak_av_mlknn, acc_all_mlknn] = mymlknn(aoos(Xtr,...
%     Database.featln, size(Xtr, 2)), Xtestorcv, cvortest, opts);

if L == 1
% KNN classifier
acc_knn_test = myknn(Xtr, Xtestorcv, Database, cvortest); % k = 5    
acc_svm_test = mysvm(Xtr, Xtestorcv, Database, cvortest);
acc_rdf_test = myrndfrst(Xtr, Xtestorcv, Database, cvortest);

% DSS_knn_test = myknn(Database.tr_data, Database.test_data, Database, cvortest); % k = 5    
% DSS_svm_test = mysvm(Database.tr_data, Database.test_data, Database, cvortest);
% DSS_rdf_test = myrndfrst(Database.tr_data, Database.test_data, Database, cvortest);
end

% dt = datestr(datetime);
% dt((datestr(dt) == ':')) = '_'; % for windows computer
% save([dt, '_test_results'], 'acc_knn_test', 'acc_svm_test','K', 'lbmd', 'mu', 'Q',...
%     'nu', 'beta', 'pctrl','mixture_n', 'seed')
else
    fprintf('no dictionary found \n')
end

toc

% run tt_esc.m % automatically run 20 shuffles
figure