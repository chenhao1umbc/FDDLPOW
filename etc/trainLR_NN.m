
% load settings
featln = 8;
lb = [1,zeros(1,5)];
for counter =1:6
    Y(featln*300*(counter-1)+1:featln*300*counter, :) = ...
    kron(ones(300*featln,1), circshift(lb, counter-1));
end
%% train LR

K = 25;
lbmd = 0.01;
mu=0.1;
nu= 10;
beta = 1;
Q= 20; % 10 for alg 1
x = [];
y = [];
for f = 1000:1004
[opts]=loadoptions(K,lbmd,mu,10,1,-1, 2000, f);
load(opts.Dictnm)
X = (Dict.W'*Dict.Z)'; % projected
x = [x; X];
y = [y; Y];
end
warning off
B = mnrfit(x, y);
save(['B_X_Y_dict1.mat'], 'B', 'X', 'Y');


lbmd = 0.01;
mu=0.1;
nu= 10;
beta = 1;
Q= 20; % 10 for alg 1
x = [];
for f = 1000:1004
[opts]=loadoptions(K,lbmd,mu,nu,1,-1, 2000, f);
load(opts.Dict2nm)
X = (Dict.W'*Dict.Z)'; % projected
x = [x; X];
end
warning off
B = mnrfit(x, y);
save(['B_X_Y_dict2.mat'], 'B', 'X', 'Y');


lbmd = 0.01;
mu=0.1;
nu= 10;
beta = 1;
Q= 20; % 10 for alg 1
x = [];
for f = 1000:1004
[opts]=loadoptions(K,lbmd,mu,nu,beta,-1, 2000, f);
load(opts.Dict3nm)
X = (Dict.W'*Dict.Z)'; % projected
x = [x; X];
end
warning off
B = mnrfit(x, y);
save(['B_X_Y_dict3.mat'], 'B', 'X', 'Y');

% % use LR
% load('dict2_k25_lmbd0.01_mu0.1_Q6_nu1000_rng1000.mat')
% load('B_X_Y_dict1.mat')
% [acc, acc_weak_av, acc_av] = lr_test(Dict, Database, Z, B, pctrl)


%% Train neural network
inp = Dict.W'*Dict.Z;
hiddenLayerSize = 10;
net = patternnet(hiddenLayerSize);
trlb = Database.tr_label;
trlb = aoos(trlb, 8, size(trlb,2));
[net,tr] = train(net,inp,trlb);
% use NN
outputs = net(z);