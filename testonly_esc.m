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

mixture_n = 2; % mixture_n classes mixture, = 1,2 (1 means non -mixture)
pctrl.db = 3; % dynamic ratio is 0 3, 5, 10, 15 db

K = 60;
lbmd = 0.025;
mu= 0.005;
Q = 0.9;% this is wq without negative
nu = 0.03;
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
Xtr = Dict.W'*Dict.Z;
% KNN classifier
acc_knn_test = myknn(Xtr, Xtestorcv, Database, cvortest); % k = 5    
acc_svm_test = mysvm(Xtr, Xtestorcv, Database, cvortest);
% dt = datestr(datetime);
% dt((datestr(dt) == ':')) = '_'; % for windows computer
% save([dt, '_test_results'], 'acc_knn_test', 'acc_svm_test','K', 'lbmd', 'mu', 'Q',...
%     'nu', 'beta', 'pctrl','mixture_n', 'seed')
end
toc
