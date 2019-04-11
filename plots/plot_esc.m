clear 
clc
close all 

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././FDDLPOW'))

%% bar chart
% zero forcing
% {'Alg.1+Zero-forcing','Alg.2+Zero-forcing','Alg.3+Zero-forcing'}
% 86.2 87.8 88.6


names = {'k-NN','Random forest','SVM','DSS+kNN', 'DSS+Random forest', 'DSS+SVM',...
    'Alg.1+k-NN','Alg.1+SVM','Alg.1+Random forest',...
    'Alg.2+k-NN','Alg.2+SVM','Alg.2+Random forest',...
    'Alg.3+k-NN','Alg.3+SVM','Alg.3+Random forest',...
    };
c =  categorical(names);
alg1_knn = [85 83 85 83 92];
alg1_rdf = [86 80 89 88 90];
alg1_svm = [90 81 90 90 92];
alg2_knn = [87 83 89 86 86];
alg2_rdf = [92 82 92 88 90];
alg2_svm = [93 84 92 88 93];
alg3_knn = [82 85 84 88 83];
alg3_rdf = [91 87 88 90 91];
alg3_svm = [91 87 87 92 92];
y = [66.7 72.7 67.5 75.8, 78.7, 78.2 ...
    mean(alg1_knn) mean(alg1_svm) mean(alg1_rdf) ...
    mean(alg2_knn) mean(alg2_svm) mean(alg2_rdf) ...
    mean(alg3_knn) mean(alg3_svm) mean(alg3_rdf) ...
    ];
bar(c, y, 0.6)
title('Non-mixture classification accuracy')
set(gcf,'color','w')

%% L=2
a1_zf = [0.6418 0.6388 0.6334 0.6081 0.5786];
a2_zf = [0.7022 0.6959 0.6884 0.6425 0.5956];
a3_zf = [0.7224 0.7161 0.7077 0.6733 0.6321]; % beta 1
% a3_zf = [0.7049 0.7011 0.6908 0.6513 0.6112]; % beta 5e-4

a1_svm = [0.5746 0.5756 0.5753 0.5689 0.5587];
a2_svm = [0.5936 0.5964 0.5916 0.5765 0.5622]; 
a3_svm = [0.5969 0.5981 0.5975 0.5767 0.5606]; % beta 5e-4


db = [0, 3, 5, 10, 15];
figure
hold on
plot(db, a1_zf, 'r--^');
plot(db, a2_zf, 'm:^');
plot(db, a3_zf, 'b-^');

plot(db, a1_svm, 'r--x');
plot(db, a2_svm, 'm:x');
plot(db, a3_svm, 'b-x');
legend('Alg.1 and zero-forcing', 'Alg.2 and zero-forcing', 'Alg.3 and zero-forcing',...
    'Alg.1 and SVM', 'Alg.2 and SVM', 'Alg.3 and SVM')
ylim([0.4, 0.75])
xticks([0, 3, 5, 10, 15])
set(gcf,'color','w');
xlabel 'Dynamic range'
ylabel('Accuracy')
title 'L = 2 classification accuracy'

%%
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


%% MDS / PCA
%{ 
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


