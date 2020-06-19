% this script is used to test LRSDL+SVM for mixture signal
close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('/extra/chenhao1/DICTOL-master/LRSDL_FDDL'))

%% load data
mixture_n = 1; % mixture_n classes mixture, = 1,2,3
SNR_INF = 2000;
pctrl.db = 20; % dynamic ratio is 0 3, 6, 10, 20 db
if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end

for f = 1:5
[Database]=load_data_new(mixture_n, SNR_INF, pctrl, f);
end
%% settings
k = 10;
k0 = 5;
lambda1 = 0.001;
lambda2 = 0.01;
lambda3 = 0.02;
C = max(label_train);

opts.k           = k;
opts.k0          = k0;
opts.show_cost   = 0;
opts.lambda1     = lambda1;
opts.lambda2     = lambda2;
opts.lambda3     = lambda3;
opts.D_range     = k*(0:C);
opts.D_range_ext = [opts.D_range k*C+k0];
opts.initmode    = 'normal';   
opts.max_iter    = 100;
opts             = initOpts(opts);
opts.verbose      = true;
opts.tol         = 1e-8;

%% Train 
[D, D0, X, X0, CoefM, coefM0, opts, rt] = LRSDL(Y_train, label_train, opts);
X1 = [X; X0];
Y_range = label_to_range(label_train);
C = max(label_train);
CoefMM0 = zeros(size(X1,1), C);
for c = 1: C 
    X1c = get_block_col(X1, c, Y_range);
    CoefMM0(:,c) = mean(X1c,2);
end    
opts.verbose = 0;
acc = [];



[acc_lrsdl, rt] = LRSDL_wrapper(Y_train, label_train, Y_test , label_test, ...
                            k, k0, lambda1, lambda2, lambda3);