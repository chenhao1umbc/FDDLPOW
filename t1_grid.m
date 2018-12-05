% table 1 

close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))

%% load settings
% do traing or do crossvalidation
do_training = 0;
cvortest = [1, 0]; % [docv, dotest] cannot be [1, 1]

mixture_n = 1; % mixture_n classes mixture, = 1,2,3 (1 means non -mixture)
pctrl.db = 0; % dynamic ratio is 0 3, 6, 10, 20 db

K = 100;
lbmd = 5e-3;
mu=1e-2;
Q = 30;% this is wq without negative
SNR = 2000;

K = [100, 200, 300];
lbmd = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 1e-4];
mu = [1, 0.1, 0.01, 0.001 0.0001];

%% load data
[Database] = load_GRID(mixture_n, SNR, pctrl);

%% training dictionary
if do_training ==1    
for ind1 = 1: length(K)
for ind2 = 1: length(lbmd)
for ind3= 1: length(mu)
    [opts]=loadoptions_grid(1,K(ind1),lbmd(ind2),mu(ind3),Q );
    % for table 1 algorithm    
    Dict = FDDLOW_table1(Database.tr_data,Database.tr_label(1,:),opts);
    if Dict.iter > 30
        save(opts.Dictnm,'Dict','opts')    
    end 
end
end
end
end

%% testing part
if sum(cvortest)
for ind1 = 1: length(K)
for ind2 = 1: length(lbmd)
for ind3= 1: length(mu)    
    [opts]=loadoptions_grid(1,K(ind1),lbmd(ind2),mu(ind3),Q );
    if exist(opts.Dictnm, 'file')        
    load(opts.Dictnm,'Dict','opts')
    Z = sparsecoding(Dict,Database,opts,mixture_n, cvortest);
    Xtestorcv = Dict.W'*Z;
    Xtr = Dict.W'*Dict.Z;
    % KNN classifier
    fprintf(opts.Dictnm)
    acc = myknn(Xtr, Xtestorcv, Database, cvortest) % k = 5        
    end
end
end
end
end

toc
