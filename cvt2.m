% table 2

close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././FDDLPOW'))
% do traing or do crossvalidation
do_training = 1;
do_result = 1;
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

% load settings
K = 25;
lbmd = 0.01;
mu=0.1;
Q=25;% this is wq without negative
nu= 10;
beta = -1;
% nu = [1e-4 5e-4 1e-3 5e-3 1e-2 5e-2 0.1 0.5];

%% testing/cv part
[Database]=load_data_new(2, SNR_INF, pctrl, 1000);
zf.acc = zeros(1,5); mf.acc = zeros(1,5);
zf.acc_weak = zeros(1,5); mf.acc_weak = zeros(1,5);
for indn = 1:length(nu)
for f = 1000:1004
if do_result ==1      
    % run doresult
    [opts]=loadoptions(K,lbmd,mu,Q,nu(indn),beta, SNR_INF, f);
    if exist(opts.Dict2nm, 'file') load(opts.Dict2nm,'Dict','opts'), else continue; end
    disp(opts.Dict2nm)
%     if exist(opts.Dictnm, 'file') load(opts.Dictnm,'Dict','opts'), else continue; end
%     disp(opts.Dictnm)
    % run prep_ZF 
    if exist('Dict')==1    Dict_mix = Dict;    end
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
    [~, zf.acc_weak(indn, f-999), zf.acc(indn, f-999)] = calc_labels(labels_pre, opts);
    [~, mf.acc_weak(indn, f-999), mf.acc(indn, f-999)] = calc_labels(labels_pre_mf, opts);
   
end
end
end
sum(zf.acc,2)/5
sum(mf.acc,2)/5
sum(zf.acc_weak,2)/5
sum(mf.acc_weak,2)/5
toc