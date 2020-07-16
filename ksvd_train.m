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


for f = 1000:1004
[Database]=load_data_new(mixture_n, SNR_INF, pctrl, f);
X = Database.tr_data;
X_cv = Database.cv_data;
label_train = Database.tr_label;
label_cv = Database.cv_label;

K = 25;
for t = [15, 18, 21]
T0 = t;
D = X(:,randi(size(X, 2), 1, K));
% parameter of the optimization procedure are chosen
param.L = T0; % not more than 10 non-zeros coefficients
param.eps=0.0; % squared norm of the residual should be less than 0.1

    for it = 1:200
    % find weights, using dictionary D 
    %If the SPAMS software (by J. Mairal et al.) is installed and available 
    % from Matlab then sparseapprox.m can be used to access the mexLasso 
    %and mexOMP functions there.             
    W = mexOMP(X,D,param);
    R = X - D*W;
    if norm(R)/norm(X) <1e-4; break; end
    for k=1:K
        I = find(W(k,:));
        Ri = R(:,I) + D(:,k)*W(k,I);
        [U,S,V] = svds(Ri,1,'L');
        D(:,k) = U;
        W(k,I) = S*V';
        R(:,I) = Ri - D(:,k)*W(k,I);
    end
    end

params = ['k_',num2str(K), 'T0_',num2str(T0), 'f_', num2str(f)];
printf(params)
save([params,'ksvd_train.mat'], 'D','W')
end
end

