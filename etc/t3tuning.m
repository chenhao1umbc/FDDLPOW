% table 3

close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././FDDLPOW'))

%% load settings
% do traing or do crossvalidation
do_training = 0;
cvortest = [1, 0]; % [docv, dotest] cannot be [1, 1]

mixture_n = 2; % mixture_n classes mixture, = 1,2,3 (1 means non -mixture)
pctrl.db = 10; % dynamic ratio is 0 3, 6, 10, 20 db
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
SNR = 2000;
beta = [7, 5, 3, 1, 0.7, 0.5, 0.3, 0.1, 0.07, 0.05, 0.03, 0.02, 0.017,...
    0.013, 0.01, 0.007, 0.005, 0.003, 0.001,7e-4, 5e-4, 3e-4, 1e-4];

%% load data
[Database] = load_ESC(mixture_n, SNR, pctrl); 
fmax =10; acc = zeros( length(beta), fmax);
for f = 1:fmax
seed = f*100;% change ramdom seed to do m-fold cv   
Database = myshuffle(Database,seed);

%% testing part
if sum(cvortest)
addpath(genpath('.././tempresult'))
for ind6 = 1:length(beta)           
    [opts]=loadoptions_ESC(3,K,lbmd,mu,Q*K,nu,beta(ind6) );
    if exist(opts.Dictnm, 'file')        
        opts.Dictnm
    load(opts.Dictnm,'Dict','opts')
    Z = sparsecoding(Dict,Database,opts,mixture_n, cvortest);
    Z = aoos(Z,Database.featln, size(Z, 2));
    Xtestorcv = Dict.W'*Z;
    Xtr = Dict.W'*Dict.Z;
    opts.n = Database.N_c; 
    opts.C = max(Database.tr_label); % 10 classe
    opts.Ncombs = max(Database.cv_mixlabel);
    opts.ln_test = size(Database.test_mixlabel, 2)/Database.featln;
    opts.equal = pctrl.equal;
    
%     % run zero-forcing
%     W = Dict.W;
%     C = opts.C;
%     N = size(Database.tr_label,2);
%     Nc = N / C;
%     featln = Database.featln;                   
%     H3 = kron(eye(C),ones(Nc, 1)/Nc); % M = Z*H3
%     M = Dict.Z*H3;
%     % zero forcing
%     H = W'*M;
%     result = pinv(H)*W'*Z;
%     [~, labels_pre] = sort(result, 1, 'descend');%     %         
%     [~, ~, acc_all] = calc_labels(labels_pre, opts);
%     acc(ind6,f) = acc_all
    
    % run svm
    [~, ~, acc_all] = mymlsvm(aoos(Xtr,...
    Database.featln, size(Xtr, 2)), Xtestorcv, cvortest, opts);
    acc(ind6,f) = acc_all; % k = 5    
    toc
    end
end
end
dt = datestr(datetime);
dt((datestr(dt) == ':')) = '_'; % for windows computer
save(['.././tempresult/m3log',dt, '_t3_results'], 'acc', 'beta', 'seed')
toc
end

toc
