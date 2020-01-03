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
SNR_INF = 2000;

% load data
mixture_n = 1; % mixture_n classes mixture, = 1,2,3

pctrl.db = 10; % dynamic ratio is 0 3, 6, 10, 20 db
if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end

%% training dictionary
% load settings
K = 25;
lbmd = 0.01;
mu=0.1;
Q=6;% this is wq without negative
nu= 0.01;
beta = -1;
nu = [1e-2];

tic
for f = 1000:1004
[Database]=load_data_new(mixture_n, SNR_INF, pctrl, f);
for indn = 1: length(nu)
    [opts]=loadoptions(K,lbmd,mu,Q,nu(indn),beta, SNR_INF, f);
    % for table 1 algorithm
    if do_training ==1
        Dict = FDDLOW_table2(Database.tr_data,Database.tr_label,opts);
        save(opts.Dict2nm,'Dict','opts')
    end
end
end
toc
