
% load settings
featln = 8;
lb = [1,zeros(1,5)];
for counter =1:6
    Y(featln*300*(counter-1)+1:featln*300*counter, :) = ...
    kron(ones(300*featln,1), circshift(lb, counter-1));
end

%% train LR_NN alg 1
K = 25;
lbmd = 0.01;
mu=0.1;
nu= 10;
beta = 1;
Q= 10; % 10 for alg 1
x = [];
y = [];

for f = 1000:1004
[opts]=loadoptions(K,lbmd,mu,Q,nu,beta,2000, f);
load(opts.Dictnm)
X = (Dict.W'*Dict.Z)'; % projected
x = [x; X];
y = [y; Y];
end
warning off
B = mnrfit(x, y);
save(['B_X_Y_dict1.mat'], 'B', 'x', 'y');
% % use LR
% load('dict2_k25_lmbd0.01_mu0.1_Q6_nu1000_rng1000.mat')
% load('B_X_Y_dict1.mat')
% [acc, acc_weak_av, acc_av] = lr_test(Dict, Database, Z, B, pctrl)


% train Neural Networks
inp = x';
hiddenLayerSize = 10;
net = patternnet(hiddenLayerSize);
[net,tr] = train(net,inp,y');
% % use NN
% outputs = net(z);
save(['NN_dict1.mat'], 'net', 'tr');


%% train LR_NN alg 2
K = 25;
lbmd = 0.1;
mu=0.001;
nu= 10;
beta = 1;
Q= 20; % 10 for alg 1
x = [];
y = [];
for f = 1000:1004
[opts]=loadoptions(K,lbmd,mu,Q,nu,beta,2000, f);
load(opts.Dict2nm)
X = (Dict.W'*Dict.Z)'; % projected
x = [x; X];
y = [y; Y];
end
warning off
B = mnrfit(x, y);
save('B_X_Y_dict2.mat', 'B', 'x', 'y');

% train Neural Networks
inp = x';
hiddenLayerSize = 10;
net = patternnet(hiddenLayerSize);
[net,tr] = train(net,inp,y');
% % use NN
% outputs = net(z);
save(['NN_dict2.mat'], 'net', 'tr');



%% train LR_NN alg 3
K = 25;
lbmd = 0.1;
mu=0.001;
nu= 10;
beta = 1;
Q= 20; % 10 for alg 1
x = [];
y = [];
for f = 1000:1004
[opts]=loadoptions(K,lbmd,mu,Q,nu,beta,2000, f);
load(opts.Dict3nm)
X = (Dict.W'*Dict.Z)'; % projected
x = [x; X];
y = [y; Y];
end
warning off
B = mnrfit(x, y);
save('B_X_Y_dict3.mat', 'B', 'x', 'y');

% train Neural Networks
inp = x';
hiddenLayerSize = 10;
net = patternnet(hiddenLayerSize);
[net,tr] = train(net,inp,y');
% % use NN
% outputs = net(z);
save(['NN_dict3.mat'], 'net', 'tr');




