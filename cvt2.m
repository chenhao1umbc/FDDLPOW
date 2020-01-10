close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././FDDLPOW'))
% do traing or do crossvalidation
cvortest = 1;  % 1 means cv, 0 means test
SNR_INF = 2000;

% load data
mixture_n = 2; % mixture_n classes mixture, = 1,2,3

pctrl.db = 0; % dynamic ratio is 0 3, 6, 10, 20 db
if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end
if mixture_n < 3  pctrl.if2weak = 0; end
K = 25;
lbmd = [0.1 0.01 0.001 0.001];
mu=[1 0.1 0.01 1e-3 1e-4] ;
Q= [6 10 20 25];
nu= [1e-3 1e-2 0.1 1 10 100];
beta = -1;

%% testing/cv part
[Database]=load_data_new(2, SNR_INF, pctrl, 1000);
zf.acc = zeros(length(mu), 7, length(nu),length(lbmd), length(Q),5); 
mf.acc = zf.acc; zf.acc_weak = zf.acc; mf.acc_weak = zf.acc;

for indll = 1: length(lbmd)
    
lamb_range = lbmd(indll)*5.^(-3:3);
for indm = 1:length(mu)
for indq = 1:length(Q)
for indn = 1:length(nu)
for f = 1000:1004    
    % run doresult
    [opts]=loadoptions(K,lbmd(indll),mu(indm),Q(indq),nu(indn),beta, SNR_INF, f);
    if exist(opts.Dict2nm, 'file') load(opts.Dict2nm,'Dict','opts'), else continue; end
    disp(opts.Dict2nm)
    for indl = 1:length(lamb_range)
        opts.lambda1 = lamb_range(indl);
    % run prep_ZF 
    Z = sparsecoding(Dict, Database, opts, mixture_n, cvortest);
    Z = aoos(Z,Database.featln, size(Z, 2));
    
    run calc_M
    % zero forcing
    H = W'*M;
    r_zf = pinv(H)*W'*Z;
    [~, labels_pre] = sort(r_zf, 1, 'descend');
    
    % matched filter
    r_matched = H'*W'*Z;
    [~, labels_pre_mf] = sort(r_matched, 1, 'descend');
    
    % calculate accuracy
    [~, zf.acc_weak(indm, indl, indn, indll, indq, f-999), zf.acc(indm, indl, indn, indll, indq, f-999)] = calc_labels(labels_pre, opts);
    [~, mf.acc_weak(indm, indl, indn, indll, indq, f-999), mf.acc(indm, indl, indn, indll, indq, f-999)] = calc_labels(labels_pre_mf, opts);
   
    end
end
end
save('alg2_0db_L2.mat','zf', 'mf','nu','beta','Q','lbmd')
end
end
end
save('alg2_0db_L2.mat','zf', 'mf','nu','beta','Q','lbmd')
nf = 5;
sum(zf.acc,3)/nf
sum(mf.acc,3)/nf
sum(zf.acc_weak,3)/nf
sum(mf.acc_weak,3)/nf