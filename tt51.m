% table 1 

close all
clear
clc;
tic

addpath(genpath('./core'));
addpath(genpath('./etc'));
addpath(genpath('/home/chenhao1/Matlab/FDDLOW/data'));
SNR_INF = 2000;

% do traing or do crossvalidation
do_training = 1;
do_cv = 0;
mixture_n = 1; % mixture_n classes mixture, = 1,2,3
pctrl.db = 0; % dynamic ratio is 0 3, 6, 10, 20 db
pctrl.if2weak = 0; % if 2 weak components in mixture of 3
if pctrl.db == 0     pctrl.equal = 1; else    pctrl.equal = 0; end
cvortest = 1;  % 1 means cv, 0 means test

%% training dictionary
% load settings
K = 50;
lbmd = 0.01;
mu=0.001;
nu= 1e3;
beta = 1;
Q= 30;
SNR = 20;

% K = [50, 100, 150, 200, 250];
lbmd = [ 0.001, 1e-4];
mu = [1, 0.1, 0.01, 0.001 0.0001];
% SNR = [2000, 20, 0, -5, -10, -20];
Q = [10 20 30 40 50 ];

for f = 1000:1004
[Database]=load_data_new(mixture_n, SNR_INF, pctrl, f);

for indk = 1: length(K)
for indq = 1: length(Q)
for indl = 1: length(lbmd)
for indm = 1: length(mu)   
    % for table 1 algorithm
    if do_training ==1
        [opts] = loadoptions(K(indk),lbmd(indl),mu(indm),Q(indq),nu,beta,SNR,f)
        Dict = FDDLOW_table1(Database.tr_data,Database.tr_label,opts);
        save(opts.Dictnm,'Dict','opts')
    end
end 
end
end
end
end
