close all
clear
clc;
tic

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
% load data
[Database]=load_data;

%% training dictionary
for ii = 32
    for jj = 1:10
        K=100;
        lbmd=1.5;
        mu=0.1;
        nu=1e3;
        Q= ii;
        beta = 1.7;
        SNR = 2000;
        [opts]=loadoptions(K,lbmd,mu,Q,nu,beta, SNR);
        opts.rng = jj;

        % for table 2 algorithm
        Dict_mix = FDDLOW_table2(Database.tr_data,Database.tr_label,opts);
        save(['rng', num2str(jj), opts.mixnm],'Dict_mix','opts')

        % % this part will give the logistic regression result
        % trainLR_dif(opts, Database);
    end
end