% table 3

close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././FDDLPOW'))
addpath(genpath('.././tempresult/'))

%% load settings
% do traing or do crossvalidation
do_training = 1;
cvortest = [1, 0]; % [docv, dotest] cannot be [1, 1]

mixture_n = 2; % mixture_n classes mixture, = 1,2,3 (1 means non -mixture)
pctrl.db = 0; % dynamic ratio is 0 3, 6, 10, 20 db

K = 60;
lbmd = 0.025;
mu= 0.005;
Q = 0.9;% this is wq without negative
nu = 0.03;
SNR = 2000;
% beta = 0.01;

% K = [20, 40, 60, 80, 100, 120 ];
% lbmd = [0.005, 0.01,0.04, 0.07, 0.1 0.4, 0.7, 1 ];
% mu = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 5e-4, 1e-4];
% Q = [1, 0.9, 0.75, 0.5, 0.3 ]; % prtion
% nu = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 5e-4, 1e-4];
beta = [7, 5, 3, 1, 0.7, 0.5, 0.3, 0.1, 0.07, 0.05, 0.03, 0.02, 0.017,...
    0.013, 0.01, 0.007, 0.005, 0.003, 0.001,7e-4, 5e-4, 3e-4, 1e-4];

%% load data
[Database] = load_ESC(mixture_n, SNR, pctrl);
acc_knn = zeros(length(K), length(lbmd), length(mu),length(Q), length(nu), length(beta));
acc_svm = zeros(length(K), length(lbmd), length(mu),length(Q), length(nu), length(beta));
for f = 1:1
seed = f*100;% change ramdom seed to do m-fold cv   
Database = myshuffle(Database,seed);
%% training dictionary
if do_training ==1    
for ind1 = 1:length(K)
for ind2 = 1:length(lbmd)
for ind3 = 1:length(mu)
for ind4 = 1:length(Q)
for ind5 = 1:length(nu)
for ind6 = 1:length(beta)
    [opts]=loadoptions_ESC(3,K(ind1),lbmd(ind2),mu(ind3),Q(ind4)*K(ind1), nu(ind5), beta(ind6) );
    % for table 1 algorithm    
    Dict = FDDLOW_table3(Database.tr_data,opts);
    if Dict.iter/opts.max_iter > 0.3
        sparsity=mean(sum(Dict.Z ~= 0))/opts.K
        save(['.././tempresult/',opts.Dictnm],'Dict','opts')    
    end 
end
end
end
end
end
end
end

%% testing part
if sum(cvortest)
addpath(genpath('.././tempresult'))
acc_knn = zeros(length(K), length(lbmd), length(mu),length(Q), length(nu), length(beta));
acc_svm = zeros(length(K), length(lbmd), length(mu),length(Q), length(nu), length(beta));
for ind1 = 1: length(K)
for ind2 = 1: length(lbmd)
for ind3= 1: length(mu)   
for ind4 = 1:length(Q)
for ind5 = 1:length(nu)
for ind6 = 1:length(beta)           
    [opts]=loadoptions_ESC(3,K(ind1),lbmd(ind2),mu(ind3),Q(ind4)*K(ind1), nu(ind5),beta(ind6) );
    if exist(opts.Dictnm, 'file')        
        opts.Dictnm
    load(opts.Dictnm,'Dict','opts')
    Z = sparsecoding(Dict,Database,opts,mixture_n, cvortest);
    Z = aoos(Z,Database.featln, size(Z, 2));
    Xtestorcv = Dict.W'*Z;
    Xtr = Dict.W'*Dict.Z;%*aoos(Dict.Z,Database.featln, size(Dict.Z, 2));
    % tuning based on non-miture
    acc_knn(ind1, ind2, ind3, ind4, ind5,ind6,f) = myknn(Xtr, Xtestorcv, Database, cvortest); % k = 5    
    acc_svm(ind1, ind2, ind3, ind4, ind5,ind6,f) = mysvm(Xtr, Xtestorcv, Database, cvortest);
    maxknn(f) = max(max(max(max(max(max(max(acc_knn)))))))
    maxsvm(f) = max(max(max(max(max(max(max(acc_svm)))))))
    end
end
end
% dt = datestr(datetime);
% dt((datestr(dt) == ':')) = '_'; % for windows computer
% save(['.././tempresult/m3log',dt, '_t3_results'], 'acc_knn', 'acc_svm', 'maxknn', 'maxsvm', 'seed')
end
end
end
end
end
end
meanknn = max(max(max(max(max(max(sum(acc_knn,7)/5))))));
meansvm = max(max(max(max(max(max(sum(acc_svm,7)/5))))));
% dt = datestr(datetime);
% dt((datestr(dt) == ':')) = '_'; % for windows computer
% save([dt, '_m3log_t3_results'], 'acc_knn', 'acc_svm', 'maxknn', 'maxsvm', 'K',...
%     'meansvm','meanknn', 'lbmd', 'mu', 'Q', 'nu', 'seed')
toc
