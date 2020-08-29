% run lrsdl.m in the main folder to add supporting directories

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


for i = [1e-4, 1e-4, 5e-5,1e-4, 1e-4, 1e-4, 1e-4,1e-4, 5e-5;
        5e-3, 5e-5, 5e-2, 5e-3, 5e-2, 1e-2, 5e-3, 5e-5, 5e-2; 
        5e-5, 0.01, 1e-4, 5e-2, 5e-5, 5e-3, 5e-4, 5e-2, 5e-4;]
lambda1     = i(1);
lambda2     = i(2);
lambda3     = i(3);

% f = 1000;
for f = 1000:1004
[Database]=load_data_new(mixture_n, SNR_INF, pctrl, f);

%% settings
k = 4;
k0 = 3;
% lambda1 = 1e-4;
% lambda2 = 5e-3;
% lambda3 = 5e-5;
C = 6;

opts.k           = k;
opts.k0          = k0;
opts.show_cost   = 0;
opts.D_range     = k*(0:C);
opts.D_range_ext = [opts.D_range k*C+k0];
opts.initmode    = 'other';   
opts.max_iter    = 100;
opts             = initOpts(opts);
opts.verbose      = false;
opts.tol         = 1e-8;

% for lambda1 = [0.005, 0.001, 0.0005, 0.0001, 5e-5]
%     for lambda2 = [0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 5e-5]
%         for lambda3 = [0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 5e-5]
opts.lambda1     = lambda1;
opts.lambda2     = lambda2;
opts.lambda3     = lambda3;



%% Train 
Y_train = Database.tr_data;
Y_cv = Database.cv_data;
label_train = Database.tr_label;
label_cv = Database.cv_label;
rmpath(genpath('/extra/chenhao1/FDDLOW/fddlow'))  % to avoid the two fista.m confusion
[D, D0, X, X0, CoefM, CoefM0, opts, rt] = LRSDL(Y_train, label_train, opts);
param = ['k_',num2str(k),'k0_',num2str(k0), 'l1_',num2str(lambda1), ...
    'l2_',num2str(lambda2), 'l3_',num2str(lambda3), 'f_', num2str(f)];
save([param,'lrscdl_train.mat'], 'D','D0','X', 'X0', 'CoefM','CoefM0', 'opts')

printf(param)
% acc1 = LRSDL_pred(Y_cv, D, D0, CoefM, CoefM0, opts, label_cv);
[Z, Z0] = local_sparse_coding(Y_cv, D, D0, CoefM0, lambda1, lambda2);
addpath(genpath('/extra/chenhao1/FDDLOW/fddlow'))
Z = aoos(Z,Database.featln, size(Z, 2)); 

acc_knn(f-999) = myknn(X, Z, Database, 1) % k = 5 ;

%         end
%     end
% end

end % end of f=1000:1004
sum(acc_knn)/5
end  % end of i 
