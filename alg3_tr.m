% table 3

close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././FDDLPOW'))

% do traing or do crossvalidation
cv = 1; % validation or testing
root = '.././data/';

% load data
mixture_n = 1; % mixture_n classes mixture, = 1,2,3
SNR_INF = 2000;
pctrl.db = 20; % dynamic ratio is 0 3, 6, 10, 20 db
if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end


%% training dictionary
% load settings
K = 25;
lbmd = 0.1;
mu=0.001;
Q=20;
nu= 10 ;
beta = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from table one we know that there are combinationes accuracy is above 0.99
% one is K = 100, lambda = 1e-4, mu = 1e-3, nu = 0.01
% another is K = 100, lambda = 1e-3, mu = 0.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beta = [1e-4,1e-3, 1e-2, 0.1, 1, 10];

for f = 1005:1009
[Database]=load_data_new(mixture_n, SNR_INF, pctrl, f);
tic
for indb = 1: length(beta)
% for table 3 algorithm
    [opts] = loadoptions(K,lbmd,mu,Q,nu,beta(indb),SNR_INF,f);
    if exist(opts.Dict3nm, 'file') continue; end
    disp(opts.Dict3nm)
    Dict = FDDLOW_table3(Database.tr_data,opts);
    toc
    save([root, opts.Dict3nm],'Dict','opts')

end
end

toc