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

mixture_n = 2; % mixture_n classes mixture, = 1,2 (1 means non -mixture)
pctrl.db = 3; % dynamic ratio is 0 3, 5, 10, 15 db
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
alg_n = 3; % algorithm number

%% load data
[Database] = load_ESC(mixture_n, SNR, pctrl);
%% cross-val part
addpath(genpath('.././tempresult'))
[opts]=loadoptions_ESC(alg_n ,K,lbmd,mu,Q*K,nu, beta);    
if exist(opts.Dictnm, 'file')
load(opts.Dictnm,'Dict','opts')
Z = sparsecoding(Dict,Database,opts,mixture_n, cvortest);
Z = aoos(Z,Database.featln, size(Z, 2));
Xtestorcv = Dict.W'*Z;
Xtr = Dict.W'*Dict.Z;

% run do_esc_zf % output acc_all acc_weak acc_weak_av
[acc_weak, acc_weak_av, acc_all] =mymlknn(Xtr, Xtestorcv, Database, cvortest);

if mixture_n == 1
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
end
toc
figure