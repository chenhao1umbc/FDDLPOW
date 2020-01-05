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
cvortest = 1;  % 1 means cv, 0 means test
SNR_INF = 2000;

% load data
mixture_n = 2; % mixture_n classes mixture, = 1,2,3

pctrl.db = 10; % dynamic ratio is 0 3, 6, 10, 20 db
if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end
if mixture_n < 3  pctrl.if2weak = 0; end
%% training dictionary
% load settings
K = 25;
lbmd = 0.001;
mu=[ 1 0.1 0.01 0.001 1e-4] ;
Q= [25 20 10 6];
nu= 0.01;
beta = -1;
%nu= [0.005, 0.01, 0.05 0.1 0.5];

% K = [50, 100, 150, 200, 250];
% lbmd = [0.1, 0.01, 0.001, 1e-4];
% mu = [1, 0.1, 0.01, 0.001 0.0001];
% SNR = [2000, 20, 0, -5, -10, -20];
% Q = [6 10 20 30 50 75 100];
for f = 1000:1004
[Database]=load_data_new(mixture_n, SNR_INF, pctrl, f);
tic
for indq = 1: length(Q)
for ind2 = 1: length(lbmd)
for ind3 = 1: length(mu)   
    % for table 1 algorithm
    if do_training ==1
        [opts] = loadoptions(K,lbmd(ind2),mu(ind3),Q(indq),nu,beta,SNR_INF,f);
        if exist(opts.Dict2nm, 'file') continue; end
        disp(opts.Dict2nm)
        Dict = FDDLOW_table2(Database.tr_data,Database.tr_label,opts);
        toc
        save(opts.Dict2nm,'Dict','opts')
    end
end 
end
end
end
toc