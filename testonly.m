
close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././FDDLPOW'))

addpath(genpath('D:\Stored_Data\data'))
SNR_INF = 2000;
cvortest = 0;

%% testing part
% load settings
K = 25;
lbmd = 0.001;
mu=0.1;
nu= 10;
beta = 1;
Q= 20; % 10 for alg 1

% %% 
% mixture_n = 1; % mixture_n classes mixture, = 1,2,3
% pctrl.db = 0; % dynamic ratio is 0 3, 6, 10, 20 db
% pctrl.if2weak = 0; % if 2 weak components in mixture of 3
% if pctrl.db == 0     pctrl.equal = 1; else    pctrl.equal = 0; end
% Database = load_data_new(mixture_n, SNR_INF, pctrl, 1000);   
% 
% r_knn = zeros(3,5,5); % 3 algs; L=1,2,3; 5 dynamic ratio; 5 folds
% SNR = [20, 0, -10, -20, -30];
% for ind_snr = 1:length(SNR)
%     db = Database;
%     db.test_data = awgn(db.test_data, SNR(ind_snr), 'measured');
%     
% for f = 1000:1004
% for alg = 1:3
%     
%     if alg ==1 Q=10; else Q=20; end
%     [opts]=loadoptions(K,lbmd,mu,Q,nu,beta, SNR_INF, f);
%     if alg == 1 load(opts.Dictnm,'Dict','opts');disp(opts.Dictnm); end
%     if alg == 2 load(opts.Dict2nm,'Dict','opts');disp(opts.Dict2nm); end
%     if alg == 3 load(opts.Dict3nm,'Dict','opts');disp(opts.Dict3nm); end
%          
%     % run prep_ZF 
%     Z = sparsecoding(Dict, db, opts, mixture_n, 0);
%     Z = aoos(Z,db.featln, size(Z, 2));
%     
%     Xtestorcv = Dict.W'*Z;
%     Xtr = Dict.W'*Dict.Z;%aoos(Dict.Z,Database.featln, size(Dict.Z, 2));
%     % KNN classifier
%     acc_knn = myknn(Xtr, Xtestorcv, db, 0); % k = 5 ;
%     r_knn(alg, ind_snr, f-999) = acc_knn;   
% end
% end
% end
% save('test_L=1.mat','r_knn');

%% 
mixture_n = 2; 

r_zf = zeros(3,5,5); %3 algs; L=2; 5 dynamic ratio; 5 folds
r_mf = r_zf;
r_zf_weak = r_zf;
r_mf_weak = r_mf;
r_lr = r_zf;
r_nn = r_zf;
dynamic_ratio = [0, 3, 6, 10, 20];    
for indd = 1:5
    pctrl.db = dynamic_ratio(indd); % dynamic ratio is 0 3, 6, 10, 20 db
    pctrl.if2weak = 0; % if 2 weak components in mixture of 3
    if pctrl.db == 0     pctrl.equal = 1; else    pctrl.equal = 0; end
    Database = load_data_new(mixture_n, SNR_INF, pctrl, 1000);   
for f = 1000:1004
for alg = 1:3

    if alg ==1 Q=10; else Q=20; end
    [opts]=loadoptions(K,lbmd,mu,Q,nu,beta, SNR_INF, f);
    if alg == 1 load(opts.Dictnm,'Dict','opts');disp(opts.Dictnm); end
    if alg == 2 load(opts.Dict2nm,'Dict','opts');disp(opts.Dict2nm); end
    if alg == 3 load(opts.Dict3nm,'Dict','opts');disp(opts.Dict3nm); end
    
    run calc_M
    Z = sparsecoding(Dict, Database, opts, mixture_n, 0);
    Z = aoos(Z,Database.featln, size(Z, 2));
    % zero forcing
    H = W'*M;
    r_zeroforcing = pinv(H)*W'*Z;
    [~, labels_pre] = sort(r_zeroforcing, 1, 'descend');

    % matched filter
    r_matched = H'*W'*Z;
    [~, labels_pre_mf] = sort(r_matched, 1, 'descend');

    % calculate accuracy
    [~, r_zf_weak(alg, indd, f-999), r_zf(alg, indd, f-999)] = calc_labels(labels_pre, opts);
    [~, r_mf_weak(alg, indd, f-999), r_mf(alg, indd, f-999)] = calc_labels(labels_pre_mf, opts);

end
end   
end




