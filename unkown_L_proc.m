% process the unknown L result
clear
clc
load('L_unknown.mat')

dynamic_ratio = [0, 3, 6, 10, 20];  
SNR_INF = 2000;
cvortest = 0;
e = 2.718281828;

% make zf, mf into probalility vectors
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

% generate the true labels
L1 = kron(eye(6), ones(1,75));
L2 = 0;
L3 = 0;



