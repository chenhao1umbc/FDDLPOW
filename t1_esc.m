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

% K = 100;
% lbmd = 0.05;
% mu= 0.05;
% Q = 1;% this is wq without negative
SNR = 2000;

K = [20, 40, 60, 80, 100, 120 ];
lbmd = [0.005, 0.01,0.04, 0.07, 0.1 0.4, 0.7, 1 ];
mu = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 5e-4, 1e-4];
Q = [1, 0.9, 0.75, 0.5, 0.3 ]; % prtion
%% load data
[Database] = load_ESC(mixture_n, SNR, pctrl);

%% training dictionary
if do_training ==1    
for ind1 = 1: length(K)
for ind2 = 1: length(lbmd)
for ind3= 1: length(mu)
for ind4 = 1:length(Q)
    [opts]=loadoptions_ESC(1,K(ind1),lbmd(ind2),mu(ind3),Q(ind4)*K(ind1) );
    % for table 1 algorithm    
    Dict = FDDLOW_table1(Database.tr_data,Database.tr_label,opts);    
    if Dict.iter/opts.max_iter > 0.3
        sparsity=mean(sum(Dict.Z ~= 0))/opts.K
%         dt = datestr(datetime);
%         dt((datestr(dt) == ':')) = '_'; % for windows computer
        save(['.././tempdict_downsampling/', opts.Dictnm],'Dict','opts')    
    end 
end
end
end
end
end

%% testing part
if sum(cvortest)
addpath(genpath('.././tempdict_downsampling'))
acc_knn = zeros(length(K), length(lbmd), length(mu),length(Q));
acc_svm = zeros(length(K), length(lbmd), length(mu),length(Q));
for ind1 = 1: length(K)
for ind2 = 1: length(lbmd)
for ind3= 1: length(mu)   
for ind4 = 1:length(Q)
    [opts]=loadoptions_ESC(1,K(ind1),lbmd(ind2),mu(ind3),Q(ind4)*K(ind1) );
    if exist(opts.Dictnm, 'file')        
    load(opts.Dictnm,'Dict','opts')
    if Dict.iter >40
        Z = sparsecoding(Dict,Database,opts,mixture_n, cvortest);
        Z = aoos(Z,Database.featln, size(Z, 2));
        Xtestorcv = Dict.W'*Z;
        Xtr = Dict.W'*Dict.Z;%aoos(Dict.Z,Database.featln, size(Dict.Z, 2));
        % KNN classifier
        acc_knn(ind1, ind2, ind3, ind4) = myknn(Xtr, Xtestorcv, Database, cvortest); % k = 5    
        acc_svm(ind1, ind2, ind3, ind4) = mysvm(Xtr, Xtestorcv, Database, cvortest);
        maxknn = max(max(max(max(acc_knn))))
        maxsvm = max(max(max(max(acc_svm))))
    end
    end
end
end
dt = datestr(datetime);
dt((datestr(dt) == ':')) = '_'; % for windows computer
save([dt, 't1_results'], 'acc_knn', 'acc_svm', 'maxknn', 'maxsvm')
end
end
end
dt = datestr(datetime);
dt((datestr(dt) == ':')) = '_'; % for windows computer
save([dt, 't1_results'], 'acc_knn', 'acc_svm', 'maxknn', 'maxsvm')
toc
