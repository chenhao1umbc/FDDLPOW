clear 
clc

addpath(genpath('./fddlow'))
addpath(genpath('./data'))

% load('DDLMDmix4_k100_lmbd1.5_mu0.1_Q16_nu1000iter_100.mat') %FDDLO
% load('FDDLOW_mix_k100_lmbd1.5_mu0.1_Q16_nu1000_beta0.1.mat')
% load('FDDLOW_mix_k100_lmbd1.5_mu0.1_Q16_nu1000_beta1.mat')
load('FDDLOW_mix_k100_lmbd1.5_mu0.1_Q16_nu1000_beta5.mat')
% load('FDDLOW_mix_k100_lmbd1.5_mu0.1_Q16_nu1000_beta10.mat')
mixture_n = 2; % mixture_n classes mixture
[Database]=load_data(mixture_n);% the equal power mixture, 400 samples per combination
if exist('Dict')==1
    Dict_mix = Dict;
end 
Z = sparsecoding_mix_test(Dict_mix, Database, opts);
W = Dict_mix.W;
C = max(Database.tr_label);
N = size(Database.tr_label,2);
Nc = N / C;
H3 = kron(eye(C),ones(Nc, 1)/Nc); % M = Z*H3
M = Dict_mix.Z*H3;
% zero forcing
H = W'*M;
result = pinv(H)*W'*aoos(Z,4,size(Z, 2));

% calculate accuracy
% because it is equal power, testing number is 
acc = 0;
comb = combnk(1:C, mixture_n);
ntestsample = size(aoos(Z,4,size(Z, 2)), 2)/mixture_n;
ncomb = size(combnk(1:C, mixture_n), 1);
nsample_percomb = ntestsample/ncomb;
[~, ind] = sort(result, 1, 'descend');
test_label = ind(1:mixture_n, :);
for i = 1: ncomb
    for ii = 1: mixture_n
        predit_labels = test_label(:, 1+nsample_percomb*(i-1):nsample_percomb*i);
        acc = acc + length(find(predit_labels == comb(i, ii)));
    end
end
acc = acc /ntestsample/mixture_n
