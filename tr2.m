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
root = '.././data/';

% load data
mixture_n = 1; % mixture_n classes mixture, = 1,2,3

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
lbmd = 0.0001;
mu = 0.1;
Q = 25;
nu= 0.01;
beta = -1;

% K = [50, 100, 150, 200, 250];
% lbmd = [0.1, 0.01, 0.001, 1e-4];
% SNR = [2000, 20, 0, -5, -10, -20];
% Q = [6 10 20 30 50 75 100];
nu = [1e-3 0.01 0.1 1 10];

for f = 1000:1004
[Database]=load_data_new(mixture_n, SNR_INF, pctrl, f);
tic
for indm = 1: length(mu)
for indn = 1: length(nu)   
    % for table 1 algorithm
    if do_training ==1
        [opts] = loadoptions(K,lbmd,mu(indm),Q,nu(indn),beta,SNR_INF,f);
        if exist(opts.Dict2nm, 'file') continue; end
        disp(opts.Dict2nm)
        Dict = FDDLOW_table2(Database.tr_data,Database.tr_label,opts);
        toc
        save([root, opts.Dict2nm],'Dict','opts')
    end
end 
end
end
toc
run tr4.m