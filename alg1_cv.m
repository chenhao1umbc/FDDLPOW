% table 1 

close all
clear
clc;
tic

addpath(genpath('./core'));
addpath(genpath('./etc'));
addpath(genpath('/home/chenhao1/Matlab/FDDLOW/data'));
SNR_INF = 2000;

% do traing or do crossvalidation
mixture_n = 1; % mixture_n classes mixture, = 1,2,3
pctrl.db = 0; % dynamic ratio is 0 3, 6, 10, 20 db
pctrl.if2weak = 0; % if 2 weak components in mixture of 3
if pctrl.db == 0     pctrl.equal = 1; else    pctrl.equal = 0; end
cvortest = 1;  % 1 means cv, 0 means test

%% training dictionary
% load settings
K = 25;
lbmd = 0.01;
mu=0.1;
nu= 1e3;
beta = 1;
Q= 10;
SNR = 20;

% K = [50, 100, 150, 200, 250];
% lbmd = [0.1, 0.01, 0.001, 1e-4];
mu = [ 0.01, 0.001 0.0001];
% SNR = [2000, 20, 0, -5, -10, -20];
% Q = [6 10 20 30 50 75 100];

%% CV/testing part
SNR = -20;
r = zeros(length(mu),length(lbmd), length(Q), length(K), 5); % Q, lambda, folds 
r_zf = r; r_mf = r;
for f = 1000:1009
% Database = load_data_new(mixture_n, SNR, pctrl, f);

[Database]=load_data_new(mixture_n, SNR_INF, pctrl, f);
Database.cv_data = awgn(Database.cv_data, SNR, 'measured');

for indk = 1: length(K)
for indq = 1: length(Q) 
for indl = 1: length(lbmd)
for indm = 1: length(mu)

    [opts]=loadoptions(K(indk),lbmd(indl),mu(indm),Q(indq),-1,-1, SNR, f);    
    if exist(opts.Dictnm, 'file') load(opts.Dictnm,'Dict','opts'), else continue; end
    disp(opts.Dictnm)
    % run prep_ZF 
    Z = sparsecoding(Dict, Database, opts, mixture_n, cvortest);
    Z = aoos(Z,Database.featln, size(Z, 2));
    
    % KNN classifier
    if Database.N_c == 1        
        Xtestorcv = Dict.W'*Z;
        Xtr = Dict.W'*Dict.Z;%aoos(Dict.Z,Database.featln, size(Dict.Z, 2));
        % KNN classifier
        acc_knn = myknn(Xtr, Xtestorcv, Database, cvortest) % k = 5 ;
        r(indm, indl, indq, indk, f-999) = acc_knn;
    end

    run calc_M
    % zero forcing
    H = W'*M;
    r_zf = pinv(H)*W'*Z;
    [~, labels_pre] = sort(r_zf, 1, 'descend');
    
    % matched filter
    r_matched = H'*W'*Z;
    [~, labels_pre_mf] = sort(r_zf, 1, 'descend');
    
    % calculate accuracy
    if Database.N_c == 1
        acc_ZF = myZFMF(labels_pre, Database, cvortest)
        acc_MF = myZFMF(labels_pre_mf, Database, cvortest)
        r_zf(indm, indl, indq, indk, f-999) = acc_ZF;
        r_mf(indm, indl, indq, indk, f-999) = acc_MF;
    else
        [acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts);
        [t, tt, ttt] = calc_labels(labels_pre_mf, opts);
    end
    
end 
end
end
end
end


toc

figure
