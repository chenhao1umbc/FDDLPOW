
% table 1 

close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
% addpath(genpath('D:\Stored_Data\data'))
addpath(genpath('.././FDDLPOW'))
SNR_INF = 2000;

% do traing or do crossvalidation
do_training = 0;
do_cv = 1;
mixture_n = 1; % mixture_n classes mixture, = 1,2,3
pctrl.db = 0; % dynamic ratio is 0 3, 6, 10, 20 db
pctrl.if2weak = 0; % if 2 weak components in mixture of 3
if pctrl.db == 0     pctrl.equal = 1; else    pctrl.equal = 0; end
cvortest = 1;  % 1 means cv, 0 means test

% load settings
K = 25;
lbmd = 0.001;
mu=0.1;
nu= 10;
beta = 1;
Q= 20;

r_knn = zeros(3,3,5,5); % 3 algs; L=1,2,3; 5 dynamic ratio; 5 folds
r_zf = r_knn; r_mf = r_knn;
dynamic_ratio = [0, 3, 6, 10, 20];


%% CV/testing part
for f = 1000:1004
for mixture_n = 1:3 
for alg = 1:3
    pctrl.db = 0; % dynamic ratio is 0 3, 6, 10, 20 db
    if pctrl.db == 0     pctrl.equal = 1; else    pctrl.equal = 0; end

    Database = load_data_new(mixture_n, SNR_INF, pctrl, f);    
    [opts]=loadoptions(K,lbmd,mu,Q,nu,beta, SNR_INF, f);
    load(opts.Dictnm,'Dict','opts')
    disp(opts.Dictnm)
         
    % run prep_ZF 
    Z = sparsecoding(Dict, Database, opts, mixture_n, cvortest);
    Z = aoos(Z,Database.featln, size(Z, 2));
    
    if mixture_n ==1
        Xtestorcv = Dict.W'*Z;
        Xtr = Dict.W'*Dict.Z;%aoos(Dict.Z,Database.featln, size(Dict.Z, 2));
        % KNN classifier
        acc_knn = myknn(Xtr, Xtestorcv, Database, cvortest); % k = 5 ;
        r_knn(alg,mixture_n, 1 ,f) = acc_knn;       
    else        
        for indd = 1:5
        pctrl.db = dynamic_ratio(indd); % dynamic ratio is 0 3, 6, 10, 20 db
        pctrl.if2weak = 0; % if 2 weak components in mixture of 3
        if pctrl.db == 0     pctrl.equal = 1; else    pctrl.equal = 0; end

        run calc_M
        % zero forcing
        H = W'*M;
        r_zeroforcing = pinv(H)*W'*Z;
        [~, labels_pre] = sort(r_zeroforcing, 1, 'descend');

        % matched filter
        r_matched = H'*W'*Z;
        [~, labels_pre_mf] = sort(r_matched, 1, 'descend');

        % calculate accuracy
        [acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts);
        [t, tt, ttt] = calc_labels(labels_pre_mf, opts);
        end
    end
end
end
end


figure
