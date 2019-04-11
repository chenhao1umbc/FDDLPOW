% table 3

close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././FDDLPOW'))

cv = 0; % validation or testing
if2weak = 1; % two weak components in the mixture
Alg_n = 2;
for mixn = 3%[2, 3]
for id = 3%[0, 3, 6, 10, 20]
% load data
mixture_n = mixn; % mixture_n classes mixture, = 1,2,3
SNR = 2000;
pctrl.db = id; % dynamic ratio is 0 3, 6, 10, 20 db
pctrl.if2weak = if2weak; % only works for mixture_n == 3
if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end
% the equal power mixture, 400 samples per combination
[Database]=load_data_new(mixture_n, SNR, pctrl);

%% load settings
K = 100;
lbmd = 1e-4;
mu=1e-3;
Q=16;% this is wq without negative
SNR = 2000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from table one we know that there are combinationes accuracy is above 0.99
% one is K = 100, lambda = 1e-4, mu = 1e-3, nu = 0.01
% another is K = 100, lambda = 1e-3, mu = 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu= 0.01 ;
if Alg_n == 2
    beta = -1;
else
    beta = 0.0001; % beta = -1, for alg2
end

[opts]=loadoptions(K,lbmd,mu,Q,nu,beta, SNR);
if SNR == 2000   
if Alg_n ==1
    load(opts.Dictnm)
else
    load(opts.mixnm) 
end
else
nm = ['SNR', num2str(SNR), opts.mixnm];
load(nm)
end
                           
%% run prep_ZF 
if exist('Dict')==1
    Dict_mix = Dict;
end
if cv == 1
    Z = sparsecoding_mix_cv(Dict_mix, Database, opts); %%%%% cv or test **************
else
    Z = sparsecoding_mix_test(Dict_mix, Database, opts);
end
W = Dict_mix.W;
C = max(Database.tr_label);
N = size(Database.tr_label,2);
Nc = N / C;
opts.C = C; % 6 classes
featln = Database.featln;
opts.n = Database.N_c;                
H3 = kron(eye(C),ones(Nc, 1)/Nc); % M = Z*H3
M = Dict_mix.Z*H3;
% zero forcing
H = W'*M;
result = pinv(H)*W'*aoos(Z,featln,size(Z, 2));
[~, labels_pre] = sort(result, 1, 'descend');

opts.Ncombs = max(Database.cv_mixlabel);
if cv == 1
    N_t = size(Database.cv_mixlabel, 2); %%%%% cv or test **************
else 
    N_t = size(Database.test_mixlabel, 2); %%%%% cv or test **************
end
opts.ln_test = N_t/featln;
opts.equal = pctrl.equal;

if if2weak ==1
[acc_weak, acc_weak_av, acc_all] = calc_labels2w(labels_pre, opts);  
else
[acc_strong, acc_weak_av, acc_all] = calc_labels(labels_pre, opts);  
end
acc_all
acc_weak_av
end
end
toc
