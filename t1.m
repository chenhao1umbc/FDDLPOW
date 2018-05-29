% table 1

close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))

do_training = 0;
do_cv = 1;

% load data
mixture_n = 2; % mixture_n classes mixture
SNR = 2000;
pctrl.equal = 0; % 1 means eqaul power, 0 non-equal
pctrl.db = 3; % dynamic ratio is 3, 6, 10, 20, 40db

% the equal power mixture, 400 samples per combination
[Database]=load_data_new(mixture_n, SNR, pctrl);
% [Database]=load_data;
% Database = load_data_spectr(1);

%% training dictionary
% load settings
K = 100;
lbmd = 0.00001;
mu=0.1;
nu= 1e3;
Q=16;% this is wq without negative
beta = 1;
SNR = 2000;

for K = [100, 150, 200, 250]
for lbmd = [0.001 0.005 0.01 1e-4]
for mu=[1, 0.1, 0.01, 0.001 0.0001]
[opts]=loadoptions(K,lbmd,mu,Q,nu,beta, SNR);

% for table 1 algorithm
if do_training ==1
    Dict = FDDLOW_table1(Database.tr_data,Database.tr_label,opts);
    if Dict.iter > 80
        save(opts.Dictnm,'Dict','opts')
    end
end 

end
end
end

%% testing part

if do_cv ==1    
    for K = [100, 150, 200, 250]
        for lbmd = [0.001 0.005 0.01 1e-4]
            for mu=[1, 0.1, 0.01, 0.001 0.0001]
                [opts]=loadoptions(K,lbmd,mu,Q,nu,beta, SNR);
                exist(opts.Dictnm)
                if exist(opts.Dictnm)
                    load(opts.Dictnm,'Dict','opts')
                end 

            end
        end
    end
end




