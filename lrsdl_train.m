% this script is used to test LRSDL+SVM for mixture signal
close all
clear
clc;
tic

addpath(genpath('/extra/chenhao1/FDDLOW/fddlow'))
addpath(genpath('/extra/chenhao1/FDDLOW/data'))
addpath(genpath('/extra/chenhao1/DICTOL-master'))

%% load data
mixture_n = 2; % mixture_n classes mixture, = 1,2,3
SNR_INF = 2000;
pctrl.db = 20; % dynamic ratio is 0 3, 6, 10, 20 db
if mixture_n < 3  pctrl.if2weak = 0; end
if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end


f = 1000;
% for f = 1000:1004
[Database]=load_data_new(mixture_n, SNR_INF, pctrl, f);

%% settings

for lambda1 = [0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]
    for lambda2 = [0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]
        for lambda3 = [0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]

k = 4;
k0 = 3;
% lambda1 = 1e-4;
% lambda2 = 1e-3;
% lambda3 = 1e-2;
C = 6;

opts.k           = k;
opts.k0          = k0;
opts.show_cost   = 0;
opts.lambda1     = lambda1;
opts.lambda2     = lambda2;
opts.lambda3     = lambda3;
opts.D_range     = k*(0:C);
opts.D_range_ext = [opts.D_range k*C+k0];
opts.initmode    = 'other';   
opts.max_iter    = 100;
opts             = initOpts(opts);
opts.verbose      = false;
opts.tol         = 1e-8;

%% Train 
Y_train = Database.tr_data;
Y_cv = Database.cv_data;
label_train = Database.tr_label;
label_cv = Database.cv_label;
rmpath(genpath('/extra/chenhao1/FDDLOW/fddlow'))  % to avoid the two fista.m confusion
[D, D0, X, X0, CoefM, CoefM0, opts, rt] = LRSDL(Y_train, label_train, opts);
X1 = [X; X0];
Y_range = label_to_range(label_train);
C = max(label_train);
CoefMM0 = zeros(size(X1,1), C);
for c = 1: C 
    X1c = get_block_col(X1, c, Y_range);
    CoefMM0(:,c) = mean(X1c,2);
end    

param = ['k_',num2str(k),'k0_',num2str(k0), 'l1_',num2str(lambda1), ...
    'l2_',num2str(lambda2), 'l3_',num2str(lambda3), 'f_', num2str(f)];
save([param,'lrscdl_train.mat'], 'D','D0','X', 'X0', 'CoefM','CoefM0', 'opts')

printf(param)
acc1 = LRSDL_pred(Y_cv, D, D0, CoefM, CoefM0, opts, label_cv);
addpath(genpath('/extra/chenhao1/FDDLOW/fddlow'))
        end
    end
end
% end
















 
