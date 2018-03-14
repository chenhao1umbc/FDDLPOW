clear 
clc

addpath(genpath('./fddlow'))
addpath(genpath('./data'))

% load data
% load('DDLMDmix4_k100_lmbd1.5_mu0.1_Q16_nu1000iter_100.mat')
load('FDDLOW_mix_k100_lmbd1.5_mu0.1_Q16_nu1000_beta0.1.mat')
% load('FDDLOW_mix_k100_lmbd1.5_mu0.1_Q16_nu1000_beta1.mat')
% load('FDDLOW_mix_k100_lmbd1.5_mu0.1_Q16_nu1000_beta10.mat')
if exist('Dict')==1
    Dict_mix = Dict; % if loading FDDL
end 
D = Dict_mix.D;
W = Dict_mix.W;
Z = Dict_mix.Z;
U = Dict_mix.U;
V = Dict_mix.V;
Delta = Dict_mix.Delta;

fln = 4; % feature length
mixture_n = 2; % mixture_n classes mixture
[Database]=load_data(mixture_n, 2000, 1);
C = max(Database.tr_label); % how many classes

Z = sparsecoding_mix_test(Dict_mix, Database, opts);
% [acc, acc_av] = lr_test_new(Dict_mix, Database, Z, B);

%% Cross-validation using KNN
%{
% sparse coding for mixture signals for L=2
Z_cvmix=sparsecoding_mix(Dict_mix,Database,opts); % testing is mixture
for ii=1:5
opts.k_neighbors=ii; % how many neighors
[~, acc_avmix(ii)]=kNN_mix(Dict_mix,Database,Z_cvmix,opts);
end

% testing
[~,opts.k_neighbors]=max(acc_avmix)
[acc, acc_avtest]=kNN_mix_test(Dict_mix,Database,Z,opts)
[~,error_table ~]=kNN_mix_test_table(Dict_mix,Database,Z,opts)
%}

%% calculate fisher/orghogonal term value
%{   
sparsity = mean(sum(Z ~= 0))
X = Database.tr_data;
N = size(Z,2);
Nc = N / C;
H1 = kron(eye(C),ones(Nc)/Nc);
H2 = ones(N)/N;
H3 = kron(eye(C),ones(Nc, 1)/Nc); % M = Z*H3
M1 = eye(N) - H1;
M2 = H1 - H2;
WtZ = W'*Z;
WtZM1 = WtZ*M1;
fWZ = norm(WtZM1, 'fro')^2 - norm(WtZ*M2, 'fro')^2 +norm(WtZ,'fro')^2; % fisher term
gWZDelta = norm(WtZ*H3 -U*Delta, 'fro')^2; % orthogonal term
cWZ = norm(WtZM1' - V, 'fro')^2; % whitening term
Loss=norm(X-D*Z,'fro')^2+opts.lambda1*sum(abs(Z(:)))+opts.mu*fWZ+opts.nu*gWZDelta + opts.beta*cWZ;

% % zero forcing result
% t = W'*m;
% svd(t) % eigen values
% (W'*m)'*W'*m;
% H = W'*m;
% result = pinv(H)*W'*Z;
%}

%% MDS / PCA
%{ 1
Cn = size(Dict_mix.Z, 2)/C; % how many per classes/combinations
temp = aoos(Dict_mix.Z(:,1:end), fln, 1200*C);
% temp = Dict_mix.Z(:,1:f:1200*n);
X0 = W'*temp;
cnt = Cn/fln; % feature length is 4
xbar = (W'*Dict_mix.Z-mean(W'*Dict_mix.Z,2));
[u,s,v] = svd(xbar*xbar');
Vq= v(:,1:3);

% load vq_40MHz_f1
% PCA = vq'*X0;
PCA =Vq'*X0;
symbolpool = {'*', 'o', 'h', 's', 'd', '^', 'p'};
c = combnk(1:C,mixture_n);
nn = size(Z,2)/size(c,1)/mixture_n; % how many samples per combo
for jj = 1:1%size(c,1)
    if jj ~=0
    h = figure(1000);
    plot3(PCA(1,1:cnt),PCA(2,1:cnt),PCA(3,1:cnt), symbolpool{1})
    grid on
    grid minor
    hold on
    for ii = 2:C        
        plot3(PCA(1,1+(ii-1)*cnt:ii*cnt),PCA(2,1+(ii-1)*cnt:ii*cnt),PCA(3,1+(ii-1)*cnt:ii*cnt),symbolpool{ii})       
    end    
    temp = aoos(Z(:,1+ nn*(jj-1):nn*jj), fln, nn);
    X = Vq'*W'*temp;    
    plot3(X(1,1:end),X(2,1:end),X(3,1:end),'x')
    legend ('ble -1','bt -2','fhss1 -3','fhss2 -4','wifi1 -5','wifi2 -6', num2str(c(jj,:)))    
%     % savefig(h,[num2str(c(jj,:)),'.fig'])    
    end
end
% close all
%}

%% shouw projected data
%{ 
% show plots

for ii = 1:C
    mc = mean(Dict_mix.Z(:,1+end/C*(ii-1):end/C*ii),2);
    temp = aoos(Dict_mix.Z, fln, size(Dict_mix.Z,2));
    proj_tr1{ii} = mc'*W*W'*temp/norm(W'*mc,2);
    figure(100);
    subplot(3,floor(C/3),ii);
    plot(proj_tr1{1, ii},'x')
    title(['projected training data--class',num2str(ii),' m\_tilde Z'])
    
    t = proj_tr1{ii};
    varProj(ii) = var(t(t>0.5)); 
    
    WtZ_c = W'*Dict_mix.Z(:,1+end/C*(ii-1):end/C*ii);
    sigma_c2{ii} = (WtZ_c - W'*mc)*(WtZ_c - W'*mc)'/(size(Dict_mix.Z,2)/C-1);
    all_sigma_c2(ii) = trace(sigma_c2{ii});
    figure(30)
    hold on
    plot(diag(sigma_c2{ii}))
    xlabel('dimension #')
    ylabel('variance of each dimension')    
    legend('ble','bt','fhss1','fhss2','wifi1','wifi2')

%     figure(31);
%     subplot(3,floor(n/3),ii);
%     t = proj_tr{1, ii};
%     histfit(t(t>0),10)
    
%     temp = aoos(Z(:,1:end/2), fln, size(Z,2)/2);
%     proj2{ii} = mc'*W*W'*temp/norm(W'*mc,2);
%     figure(40);
%     subplot(3,floor(n/3),ii);
%     plot(proj2{1, ii},'x')
%     grid minor
%     title(['projected miture testing data--class',num2str(ii),' m\_tildeW^TZ'])
%     xlabel('samples number')
     
%     figure(41);
%     subplot(3,floor(n/3),ii);
%     t = proj2{1, ii};
%     histfit(t(t>0),10)    
end

figure
plot(varProj)
xlabel('class#')
ylabel('projected variance')

figure
plot(all_sigma_c2)
xlabel('class#')
ylabel('sum of covariance diagnal')
%}