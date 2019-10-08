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
cv = 1; % validation or testing

for uuu = 2000%[20, 0, -5, -10, -20, -30]
% load data
mixture_n = 1; % mixture_n classes mixture, = 1,2,3
SNR = uuu; %SNR 2000, 20, 0, -5, -10, -20, -30
pctrl.db = 10; % dynamic ratio is 0 3, 6, 10, 20 db
if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end

% the equal power mixture, 400 samples per combination
[Database]=load_data_new(mixture_n, SNR, pctrl);


%% training dictionary
% load settings
K = 100;
lbmd = 1e-4;
mu=1e-3;
Q=16;% this is wq without negative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from table one we know that there are combinationes accuracy is above 0.99
% one is K = 100, lambda = 1e-4, mu = 1e-3
% another is K = 100, lambda = 1e-3, mu = 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nu= 0.01;
beta = -1;
% nu = [ 0.003 0.005 0.007 0.01 0.03 0.05 0.07 0.1];
   
for ind1 = 1: length(nu)
    [opts]=loadoptions(K,lbmd,mu,Q,nu(ind1),beta, SNR);
    % for table 1 algorithm
    if do_training ==1
        Dict_mix = FDDLOW_table2(Database.tr_data,Database.tr_label,opts);
        if Dict_mix.iter > 30
            save(['SNR', num2str(SNR), opts.mixnm],'Dict_mix','opts')
        end
    end
end


%% testing/cv part
if do_result ==1      
    run doresult
%     save(['SNR', num2str(SNR),'tb2_results'],'result_nu','result_nuWEEK','sparsity_nu','tr_sparsity_nu')
end


end
toc

