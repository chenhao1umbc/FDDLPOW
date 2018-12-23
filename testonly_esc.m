% table 1, 2,3 test data only no training

close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))

%% load settings
% do traing or do crossvalidation
cvortest = [0, 1]; % [docv, dotest] cannot be [1, 1]

mixture_n = 1; % mixture_n classes mixture, = 1,2,3 (1 means non -mixture)
pctrl.db = 0; % dynamic ratio is 0 3, 6, 10, 20 db

K = 40;
lbmd = 0.04;
mu= 0.005;
Q = 0.5;
nu = 0.1;
beta = 0.01;
SNR = 2000;
table_n = 1; % algorithm number

%% load data
[Database] = load_ESC(mixture_n, SNR, pctrl);

%% cross-val part
addpath(genpath('.././tempresult'))
[opts]=loadoptions_ESC(table_n ,K,lbmd,mu,Q*K,nu, beta);    
if exist(opts.Dictnm, 'file')
load(opts.Dictnm,'Dict','opts')
Z = sparsecoding(Dict,Database,opts,mixture_n, cvortest);
Z = aoos(Z,Database.featln, size(Z, 2));
Xtestorcv = Dict.W'*Z;
Xtr = Dict.W'*Dict.Z;%aoos(Dict.Z,Database.featln, size(Dict.Z, 2));
% KNN classifier
acc_knn_test = myknn(Xtr, Xtestorcv, Database, cvortest); % k = 5    
acc_svm_test = mysvm(Xtr, Xtestorcv, Database, cvortest);
dt = datestr(datetime);
dt((datestr(dt) == ':')) = '_'; % for windows computer
save([dt, '_test_results'], 'acc_knn_test', 'acc_svm_test','K', 'lbmd', 'mu', 'Q',...
    'nu', 'beta', 'pctrl','mixture_n', 'seed')
end
toc
