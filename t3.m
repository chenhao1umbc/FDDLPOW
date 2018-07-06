% table 3

close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))


% do traing or do crossvalidation
do_training = 0;
do_result = 1;
cv = 1; % validation or testing


% load data
mixture_n = 3; % mixture_n classes mixture, = 1,2,3
SNR = 2000;
pctrl.db = 20; % dynamic ratio is 0 3, 6, 10, 20 db
if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end

% the equal power mixture, 400 samples per combination
[Database]=load_data_new(mixture_n, SNR, pctrl);
% [Database]=load_data;
% Database = load_data_spectr(1);

%% training dictionary
% load settings
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
beta = 0; % or 0.14
beta = [1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1 10 100];

if do_training ==1
    for ind1 = 1: length(beta)
        [opts]=loadoptions(K,lbmd,mu,Q,nu,beta(ind1), SNR);
        % for table 1 algorithm
        Dict_mix = FDDLOW_table3(Database.tr_data,Database.tr_label,opts);
        if Dict_mix.iter > 30
            save(['SNR', num2str(SNR), opts.mixnm],'Dict_mix','opts')
        end
    end
end

%% testing part
if do_result ==1      
    run doresult
    save('tb3_results20dbNc3','result_beta','result_betaWEEK','sparsity_beta','tr_sparsity_beta')
end

 

toc