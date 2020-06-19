% this script is used to test LRSDL+SVM for mixture signal
close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('/extra/chenhao1/DICTOL-master/LRSDL_FDDL'))

%% load data
mixture_n = 2; % mixture_n classes mixture, = 1,2,3
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
C = 6;

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

acc1 = LRSDL_pred(Y_test, D, D0, CoefM, coefM0, opts, label_test);


















function acc = LRSDL_pred(Y, D, D0, CoefM, m0, opts, label_test)
nClasses = size(CoefM, 2);
k = opts.k;
k0 = opts.k0;
D_range = k*(0:nClasses);

N = size(Y,2);
acc = [];
% --------------- Sparse coding -------------------------
for lambda1 = [0.001]
    [X, X0] = local_sparse_coding(Y, D, D0, m0, lambda1, lambda2);
    % --------------- classification -------------------------
    Yhat = Y - D0*X0;
    E1 = zeros(nClasses, N);
    E2 = E1;
    for c = 1: nClasses            
        Dc = get_block_col(D, c, D_range);
        Xc = get_block_row(X, c, D_range);
        Mc = repmat(CoefM(:, c), 1, N );
        R1 = Yhat - Dc*Xc;
        R2 = X - Mc;
        E1(c,:) = sum(R1.^2);
        E2(c,:) = sum(R2.^2);
    end
    for w = [ .5]
        E = w*E1 + (1-w)*E2;
        [~, pred] = min(E);
        aaaa = double(sum(pred == label_test))/N;
        acc = [acc aaaa];
        fprintf('w: %f, lambda1 = %.4f,  acc: %f\n', w, lambda1, aaaa);
    end 
end 

    function [X, X0] = local_sparse_coding(Y, D, D0, m0, lambda1, lambda2)
        N      = size(Y,2);
        k      = size(D,2);
        k0     = size(D0,2);
        X1init = zeros(k + k0, N);
        D1     = [D D0];
        M0     = repmat(m0, 1, N);
        D1tD1  = D1'*D1;
        D1tY   = D1'*Y;
        %% cost
        function cost = calc_F(X1)
            X = X1(1: k, :);
            X0 = X1(k+1:end,:);
            cost =  0.5*normF2(Y - D1*X1) + ...
                    0.5*lambda2*normF2(X0 - M0) + ...
                    lambda1*norm1(X1);
        end 
        %% grad
        function g = grad(X1)
            X  = X1(1: k, :);
            X0 = X1(k+1:end,:);
            g  = (D1tD1*X1 - D1tY + lambda2* [zeros(k, N); X0 - M0]);
        end     
        %% ========= Main FISTA ==============================
        L             = max(eig(D1tD1)) + 2;
        opts.tol      = 1e-8;
        opts.max_iter = 300;
        [X1, ~]       = fista(@grad, X1init, L, lambda1, opts, @calc_F);  
        X             = X1(1: k, :);
        X0            = X1(k+1:end,:);
    end 
end 
