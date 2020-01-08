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

pctrl.db = 20; % dynamic ratio is 0 3, 6, 10, 20 db
if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end
if mixture_n < 3  pctrl.if2weak = 0; end
K = 25;
lbmd = 0.001;
mu=0.1;
Q= 20;
nu= 10;
beta = [1e-3 1e-2 1e-1 1 10 100];

%% testing/cv part
[Database]=load_data_new(2, SNR_INF, pctrl, 1000);
lamb_range = lbmd*5.^(-3:3);

zf.acc = zeros(length(beta),length(lamb_range),5); mf.acc = zf.acc;
zf.acc_weak = zf.acc; mf.acc_weak = zf.acc;

for indb = 1:length(beta)
for f = 1000:1004     
for indl = 1:length(lamb_range)
    % run doresult
    [opts]=loadoptions(K,lbmd,mu,Q,nu,beta(indb), SNR_INF, f);
    if exist(opts.Dict3nm, 'file') load(opts.Dict3nm,'Dict','opts'), else continue; end
    disp(opts.Dict3nm)
    opts.lambda1 = lamb_range(indl);
    % run prep_ZF 
    if exist('Dict')==1    Dict_mix = Dict;    end
    Z = sparsecoding(Dict, Database, opts, 2, cvortest);
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
    [~, zf.acc_weak(indb,indl, f-999), zf.acc(indb,indl, f-999)] = calc_labels(labels_pre, opts);
    [~, mf.acc_weak(indb,indl, f-999), mf.acc(indb,indl, f-999)] = calc_labels(labels_pre_mf, opts);
end
end
end
nf = 5;
sum(zf.acc,3)/nf
sum(mf.acc,3)/nf
sum(zf.acc_weak,3)/nf
sum(mf.acc_weak,3)/nf
save('alg3_20db_L2.mat','zf', 'mf')
toc