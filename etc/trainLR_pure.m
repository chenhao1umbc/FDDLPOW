% train logistic regression for pure signal
clear 
clc
tic
addpath(genpath('./fddlow'))
addpath(genpath('./data'))

mixture_n = 1; % mixture_n classes mixture
featln = 4;
SNR = -40;

load(['SNR', num2str(SNR), 'DDLMDmix4_k100_lmbd1.5_mu0.1_Q16_nu1000iter_100.mat']);

% loaddata-load the original power; loaddata2-equal power; loadata3-10MHz  
[Database]=load_data(mixture_n, SNR);% the equal power mixture, 400 samples per combination
Z_tr = sparsecoding_train(Dict_mix, Database, opts);

X_tr = (Dict_mix.W'*Z_tr)'; % projected

lb = [1,zeros(1,5)];
for counter =1:6
    Y_tr(featln*300*(counter-1)+1:featln*300*counter, :) = ...
    kron(ones(300*featln,1), circshift(lb, counter-1));
end

%% train logistic regression
warning off
B = mnrfit(X_tr, Y_tr);
save(['SNR', num2str(SNR), 'B_X_Y_pure.mat'], 'B', 'X_tr', 'Y_tr');


[acc, acc_av] = lr_test(Dict_mix, Database, Z, B);
toc