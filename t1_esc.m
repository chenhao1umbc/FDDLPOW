% table 1 

close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))

%% load settings
% do traing or do crossvalidation
do_training = 1;
cvortest = [1, 0]; % [docv, dotest] cannot be [1, 1]

mixture_n = 1; % mixture_n classes mixture, = 1,2,3 (1 means non -mixture)
pctrl.db = 0; % dynamic ratio is 0 3, 6, 10, 20 db

K = 100;
lbmd = 1e-16;
mu=1e-2;
Q = 100;% this is wq without negative
SNR = 2000;

% K = [5, 10, 15, 20];
% lbmd = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 5e-4, 1e-4];
% mu = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 5e-4, 1e-4];
% Q = [30, 70, 100, 150];
%% load data
[Database] = load_ESC(mixture_n, SNR, pctrl);

%% training dictionary
if do_training ==1    
for ind1 = 1: length(K)
for ind2 = 1: length(lbmd)
for ind3= 1: length(mu)
for ind4 = 1:length(Q)
    [opts]=loadoptions_ESC(1,K(ind1),lbmd(ind2),mu(ind3),Q(ind4) );
    % for table 1 algorithm    
    Dict = FDDLOW_table1(Database.tr_data,Database.tr_label(1,:),opts);
    if Dict.iter > 40
        save(opts.Dictnm,'Dict','opts')    
    end 
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
for ind4 = 1:length(Q)
    [opts]=loadoptions_grid(1,K(ind1),lbmd(ind2),mu(ind3),Q(ind4) );
    if exist(opts.Dictnm, 'file')        
        opts.Dictnm
    load(opts.Dictnm,'Dict','opts')
    Z = sparsecoding(Dict,Database,opts,mixture_n, cvortest);
    Xtestorcv = Dict.W'*Z;
    Xtr = Dict.W'*Dict.Z;
    % KNN classifier
    fprintf(opts.Dictnm)
    acc = myknn(Xtr, Xtestorcv, Database, cvortest) % k = 5    
    acc = mysvm(Xtr, Xtestorcv, Database, cvortest)
    end
end
end
end
end
end

toc
