% process the unknown L result
clear
clc
addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././FDDLPOW'))
load('L_unknown.mat')

dynamic_ratio = [0, 3, 6, 10, 20];  
SNR_INF = 2000;
cvortest = 0;
e = 2.718281828;

%% make zf, mf into probalility vectors
for mixture_n = 1:3
for indd = [1,5]   
for f = 1000:1009
for alg = 1:3
    if alg == 1 && f<1005      
        % ZF detector        
        s = e.^r_zf{mixture_n, indd, f-999, alg};
        r_zf{mixture_n, indd, f-999, alg} = s./sum(s,1);

        % matched filter
        s = e.^r_mf{mixture_n, indd, f-999, alg};
        r_mf{mixture_n, indd, f-999, alg} = s./sum(s,1);
    end    
    if alg == 2 && f<1005
        % ZF detector        
        s = e.^r_zf{mixture_n, indd, f-999, alg};
        r_zf{mixture_n, indd, f-999, alg} = s./sum(s,1);

        % matched filter
        s = e.^r_mf{mixture_n, indd, f-999, alg};
        r_mf{mixture_n, indd, f-999, alg} = s./sum(s,1);
    end    
    if alg == 3 && f>1004
        % ZF detector        
        s = e.^r_zf{mixture_n, indd, f-999, alg};
        r_zf{mixture_n, indd, f-999, alg} = s./sum(s,1);

        % matched filter
        s = e.^r_mf{mixture_n, indd, f-999, alg};
        r_mf{mixture_n, indd, f-999, alg} = s./sum(s,1);
    end
end
end
end
end

% get the averaged result over the f/fold
ar_zf = cell(3,2,3); % mixture_n, 0dB 20dB, alg
ar_mf = ar_zf; ar_lr = ar_zf;  ar_nn = ar_zf;

ar_zf = avv(r_zf, ar_zf);
ar_mf = avv(r_mf, ar_mf);
ar_nn = avv(r_nn, ar_nn);
ar_lr = avv(r_lr, ar_lr);

%% generate the true labels
L1 = kron(eye(6), ones(1,75));

l2 = zeros(6, 1500);
a = combnk(1:6, 2);
nc = size(a, 1);
for i = 1: nc
   ind = a(i, :);
   l2(ind(1), ((i-1)*100+1) :i*100) = ones(1, 100);   
   l2(ind(2), ((i-1)*100+1) :i*100) = ones(1, 100); 
end
L2 = [l2, l2];


l3 = zeros(6, 2000);
a = combnk(1:6, 3);
nc = size(a, 1);
for i = 1: nc
   ind = a(i, :);
   l3(ind(1), ((i-1)*100+1) :i*100) = ones(1, 100);   
   l3(ind(2), ((i-1)*100+1) :i*100) = ones(1, 100); 
   l3(ind(3), ((i-1)*100+1) :i*100) = ones(1, 100);
end
L3 = [l3, l3, l3];

%% calculate the acc, ROC and acc vs threhold

% plot acc
% L = 2, zf, 0db
res = zeros(3, 100);
L = 2;
for alg = 1:3
    for thr = 0.01: 0.01: 1
        r = ar_zf{L, 1, alg};
        r(r>thr) = 1;
        r(r<=thr) = 0;
        rr = (L2 - r);
        res(alg, round(thr*100)) = sum(rr ==0, 'all' ) /numel(rr);
    end
end
figure(100)
plot(res')
legend

% L = 2, zf, 20db
res = zeros(3, 100);
for alg = 1:3
    for thr = 0.01: 0.01: 1
        r = ar_zf{L, 2, alg};
        r(r>thr) = 1;
        r(r<=thr) = 0;
        rr = (L2 - r);
        res(alg, round(thr*100)) = sum(rr ==0, 'all' ) /numel(rr);
    end
end
figure(101)
plot(res')
legend


% plot ROC
% L = 3, zf,0db
for alg = 1:3
r = ar_zf{2, 1, alg};
end

% L = 3, mf,0db

% L = 3, lr,0db

% L = 3, nn,0db









