% process the unknown L result
clear
clc
addpath(genpath('../../fddlow'))
addpath(genpath('../../data'))
addpath(genpath('../../FDDLPOW'))


dynamic_ratio = [0, 3, 6, 10, 20];  
SNR_INF = 2000;
cvortest = 0;
e = 2.718281828;

%% make zf, mf into probalility vectors
table = cell(3, 2, 5, 3, 4); %3 by 2 by 5 by3 by mixture_n, [0, 20], f, alg, zf/mf...
for alg = 1:3
    load(['L_unknown_alg', num2str(alg), '.mat'])
for mixture_n = 1:3
for indd = [1,5]   
for f = 1:5 
        % ZF detector        
        s = e.^r_zf{mixture_n, indd, f};
        r_zf{mixture_n, indd, f} = s./sum(s,1);

        % matched filter
        s = e.^r_mf{mixture_n, indd, f};
        r_mf{mixture_n, indd, f} = s./sum(s,1);
        
        table{mixture_n, round((indd+3)/4), f, alg, 1} = r_zf{mixture_n, indd, f};
        table{mixture_n, round((indd+3)/4), f, alg, 2} = r_mf{mixture_n, indd, f};
        table{mixture_n, round((indd+3)/4), f, alg, 3} = r_lr{mixture_n, indd, f};
        table{mixture_n, round((indd+3)/4), f, alg, 4} = r_nn{mixture_n, indd, f};
end
end
end
end

% logistic regression may give Nan which should be 1
for mixture_n = 1:3
    for indd = [1,2]   
for f = 1:5 
    for alg = 1:3
for icls = 1:4
    t = table{mixture_n, indd, f, alg, icls};
    t(isnan(t)) = 1;
    table{mixture_n, indd, f, alg, icls} = t; 
end
    end
end
    end
end


% get the averaged tableult over the f/fold
ar_zf= cell(3,2,3); % mixture_n, 0dB 20dB, alg
ar_mf = ar_zf; ar_lr = ar_zf;  ar_nn = ar_zf;

ar_zf = avv(table(:,:,:,:,1), ar_zf);
ar_mf = avv(table(:,:,:,:,2), ar_mf);
ar_lr = avv(table(:,:,:,:,3), ar_lr);
ar_nn = avv(table(:,:,:,:,4), ar_nn);

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


%%
% _________________________acc vs threhold_______________________________________%
% 
db = [0,20];
cls = {'ZF detector', 'MF detector', 'LR classifier', 'NN classifier'};
lgd = {'Algorithm 1', 'Algorithm 2', 'Algorithm 3'};
dr = 0.01: 0.01: 1;

for L = 1:3
    for idb = 1:2
        if L == 1 && idb ==2  continue; end
        for icls = 1:4

if icls == 1  ar = ar_zf; end
if icls == 2  ar = ar_mf; end
if icls == 3  ar = ar_lr; end
if icls == 4  ar = ar_nn; end
res = zeros(3, 100);
for alg = 1:3
    for thr = 0.01: 0.01: 1
        r = ar{L, idb, alg};
        r(r>thr) = 1;
        r(r<=thr) = 0;
        rr = (eval(['L', num2str(L)]) - r);
        res(alg, round(thr*100)) = sum(rr ==0, 'all' ) /numel(rr);
    end
end
figure(100*L + 10*idb + icls)
set(gcf,'color','w')
plot(dr,res(1,:), ':b',dr,res(2,:),'--r',dr, res(3,:),'-g','linewidth', 3, 'MarkerSize', 5);
ylim([0,1]);
xlabel('Threshold')
set(gca,'FontSize',14)
set(gcf, 'Position',  [100, 100, 700, 600])
ylabel('Classification accuracy')
legend(lgd)
if idb ==1 pw_rt = ', equal power'; else pw_rt = ', 20 dB power ratio'; end
title([' L = ', num2str(L), ',', cls{icls}, pw_rt])

        end
    end
end


% _________________________plot ROC_______________________________________%

fl = L3(:);
linewidth = 4;
for L = 1:3
    if L == 1  fl = L1(:); end
    if L == 2  fl = L2(:); end
    if L == 3  fl = L3(:); end
        
    for idb = 1:2
        if L == 1 && idb ==2  continue; end
        for icls = 1:4
            
if icls == 1  ar = ar_zf; end
if icls == 2  ar = ar_mf; end
if icls == 3  ar = ar_lr; end
if icls == 4  ar = ar_nn; end
figure(1000 + 100*L + 10*idb + icls)
for alg = 1:3
    r = ar{L, idb, alg};
    fr = r(:);
    [tpr, fpr, ~] = roc(fl', fr');
    if alg == 1 plot(fpr, tpr,':b', 'linewidth', linewidth); end
    if alg == 2 plot(fpr, tpr,'--r', 'linewidth', linewidth ); end
    if alg == 3 plot(fpr, tpr,'-g', 'linewidth', linewidth ); end
    hold on
end
set(gcf,'color','w')
xlabel('False positive rate')
set(gca,'FontSize',14)
set(gcf, 'Position',  [100, 100, 700, 600])
ylabel('True positive rate')
legend(lgd)
if idb ==1 pw_rt = ', equal power'; else pw_rt = ', 20 dB power ratio'; end
title(['ROC, L = ', num2str(L), ',', cls{icls}, pw_rt])

        end
    end
end


%% Calc AUC (Area under curve)





