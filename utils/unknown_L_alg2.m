% using sparse bayesian RVM to avoid the lamb_bar tuning
clear
clc
addpath(genpath('../../fddlow'))
addpath(genpath('../../data'))
addpath(genpath('../../FDDLPOW'))

dynamic_ratio = [0, 3, 6, 10, 20];  
SNR_INF = 2000;
cvortest = 0;

r_zf = cell(3, 2, 5);
r_mf = r_zf; r_lr = r_zf; r_nn = r_zf;

for mixture_n = 1:3
for indd = [1,5]
    pctrl.db = dynamic_ratio(indd); % dynamic ratio is 0 3, 6, 10, 20 db
    if mixture_n == 3  pctrl.if2weak = 1; else pctrl.if2weak = 0; end
    if pctrl.db == 0     pctrl.equal = 1; else    pctrl.equal = 0; end
    if mixture_n == 3 && pctrl.db == 0  pctrl.if2weak = 0; end
    Database = load_data_new(mixture_n, SNR_INF, pctrl, 1000);
    
for f = 1000:1004
for alg = 1:3    
    if alg == 2 && f<1005
        load(['dict2_k25_lmbd0.1_mu0.001_Q20_nu10_rng',num2str(f),'.mat']);
        load('B_X_Y_dict2.mat')
        load('NN_dict2.mat')
        if mixture_n ~= 1  opts.lambda1 = 0.035; else opts.lambda1 = 0.05;end
        disp(opts.Dict2nm); 
        
        run calc_M
%         Z = sparsecoding(Dict, Database, opts, mixture_n, cvortest); % regular sparse coding
        Z = sbrvm(Dict, database, mixture_n, cvortest);        
        Z = aoos(Z, Database.featln, size(Z, 2));
        wtz = W'*Z; %(W'*aoos(Z, featln, N));
        H = W'*M;
        
        % ZF detector
        r_zf{mixture_n, indd, f-999} = pinv(H)*wtz;

        % matched filter
        r_mf{mixture_n, indd, f-999} = H'*wtz;

        %neural networks
        r_nn{mixture_n, indd, f-999} = net(wtz);

        %logistic regression classifier
        pre_prob = mnrval(B, wtz');
        r_lr{mixture_n, indd, f-999} = pre_prob';
    end
    
end
end
end
end

data_structure = ['3 by 2 by 5  by 6*450 cell, mixture_n = 1:3,'...
    'dynamic_ratio = [0, 20], f = 1000:1004,'];
save('L_unknown_alg2.mat','r_zf', 'r_mf',  'r_nn','r_lr','data_structure');
figure;