% clear 
% clc
% close all 

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))
addpath(genpath('.././FDDLPOW'))


%% MDS / PCA
% { 
fln = 1; % 1 means no average over one sample
C = 6;
mixture_n = 1;
Cn = size(Dict.Z, 2)/C; % how many sampless per classes/combinations
temp = aoos(Dict.Z(:,1:end), fln, size(Dict.Z, 2));
% temp = Dict.Z(:,1:f:1200*n);
W = Dict.W;
X0 = W'*temp;
cnt = Cn/fln; % feature length is 4
xbar = (W'*Dict.Z-mean(W'*Dict.Z,2));
[u,s,v] = svd(xbar*xbar');
Vq= v(:,1:3);

PCA =Vq'*X0;
symbolpool = {'*', 'o', 'h', 's', 'd', '^', 'p', '.', '+', 'x', '>'};
c = combnk(1:C,mixture_n);
% nn = size(Z,2)/size(c,1)/mixture_n; % how many samples per combo
axislimits = [-1 1 -1 1 -1 1];
jjmax = 1;
for jj = 1:1
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
%     if plottesting == 1
%         temp = aoos(Z(:,1+ nn*(jj-1):nn*jj), fln, nn);
%         X = Vq'*W'*temp;    
%         plot3(X(1,1:end),X(2,1:end),X(3,1:end),'x')
%     end
    legend ('Chainsaw','Clock tick','Cracking fire','Crying baby',...
        'Dog','Helicopter','Rain','Rooster','Sea waves','Sneezing', num2str(c(jj,:)))    
    % savefig(h,[num2str(c(jj,:)),'.fig'])    
    set(gcf,'color','w');
end
% close all
%}


