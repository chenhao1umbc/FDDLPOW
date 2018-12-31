clear 
clc


addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))

plottesting =0; % if plottesing == 1, the plot teing samples
SNR = 0; % SNR could be 2000,20, 0, -5, -10, -20, -30
fln = 1; % feature length
mixture_n = 1; % mixture_n classes mixture, = 1,2,3
cvortest = [0, 1]; % [docv, dotest] cannot be [1, 1]
pctrl.db = 3; % dynamic ratio is 0 3, 5, 10, 15 db
K = 60;
lbmd = 0.025;
mu= 0.005;
Q = 0.9;% this is wq without negative
nu = 0.03;
beta = 0.01;
SNR = 2000;
table_n = 3; % algorithm number

%% load data
[opts]=loadoptions_ESC(table_n ,K,lbmd,mu,Q*K,nu, beta);    
load(opts.Dictnm,'Dict','opts')
if exist('Dict')==1
    Dict_mix = Dict; % if loading FDDL
end 
D = Dict_mix.D;
W = Dict_mix.W;
Z = Dict_mix.Z;
if isfield(Dict_mix, 'U')
    U = Dict_mix.U;
end
if isfield(Dict_mix, 'V')
    V = Dict_mix.V;
    Delta = Dict_mix.Delta;
end
[Database] = load_ESC(mixture_n, SNR, pctrl);
C = max(Database.tr_label); % how many classes
if plottesting == 1
 Z = sparsecoding_mix_cv(Dict_mix, Database, opts);
end


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
H_bar_i = M1(1:Nc, 1:Nc);

WtZ = W'*Z;
WtZM1 = WtZ*M1;
fWZ = norm(WtZM1, 'fro')^2 - norm(WtZ*M2, 'fro')^2 + 1.1*norm(WtZ,'fro')^2; % fisher term
if isfield(Dict_mix, 'U')
    gWZ = norm(WtZ*H3 -U, 'fro')^2; % orthogonal term
    Loss=norm(X-D*Z,'fro')^2+opts.lambda1*sum(abs(Z(:)))+opts.mu*fWZ+opts.nu*gWZ;
end

Projected_WithinClass_covariance = norm(WtZM1*WtZM1','fro')
Frobenius_norm_W = norm(W, 'fro')
%}

%% MDS / PCA
%{1
Cn = size(Dict_mix.Z, 2)/C; % how many sampless per classes/combinations
temp = aoos(Dict_mix.Z(:,1:end), fln, size(Dict_mix.Z, 2));
% temp = Dict_mix.Z(:,1:f:1200*n);
X0 = W'*temp;
cnt = Cn/fln; % feature length is 4
xbar = (W'*Dict_mix.Z-mean(W'*Dict_mix.Z,2));
[u,s,v] = svd(xbar*xbar');
Vq= v(:,1:3);

PCA =Vq'*X0;
symbolpool = {'*', 'o', 'h', 's', 'd', '^', 'p', '.', '+', 'x', '>'};
c = combnk(1:C,mixture_n);
nn = size(Z,2)/size(c,1)/mixture_n; % how many samples per combo
if plottesting == 1
    jjmax = size(c,1);
else
    jjmax = 1;
end
axislimits = [-1 1 -1 1 -1 1];
for jj = 1:jjmax
    figure
    plot3(PCA(1,1:cnt),PCA(2,1:cnt),PCA(3,1:cnt), symbolpool{1})
    grid on
    grid minor
    hold on
    for ii = 2:C        
        plot3(PCA(1,1+(ii-1)*cnt:ii*cnt),PCA(2,1+(ii-1)*cnt:ii*cnt),PCA(3,1+(ii-1)*cnt:ii*cnt),symbolpool{ii})       
    end    
%     axis(axislimits)
    %  plot the testing data
    if plottesting == 1
        temp = aoos(Z(:,1+ nn*(jj-1):nn*jj), fln, nn);
        X = Vq'*W'*temp;    
        plot3(X(1,1:end),X(2,1:end),X(3,1:end),'x')
    end
    legend ('Chainsaw','Clock tick','Cracking fire','Crying baby',...
        'Dog','Helicopter','Rain','Rooster','Sea waves','Sneezing', num2str(c(jj,:)))    
    % savefig(h,[num2str(c(jj,:)),'.fig'])    
    set(gcf,'color','w');
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

