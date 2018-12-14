function [D,Z,W,Loss,opt]=initdict(X,trlabels,opt)
% This fucntio is to initialize Dictionary
% The input is  X, the training data, a matrix M by N, N data samples
%             trlabels is the labels of training data, like[1,1,1,1,2,2,3,3,3]             
%             opt are the training options with
%                 opt.K -the number of atoms in Dictionary
%                 opt.Q -the projected dimensions
%                 opt.lambda1 -lambda1 control the sparsity level
%                 opt.mu -mu the coeffecient for fisher term
%                 opt.max_iter - the maximum iteration
%                 opt.losscalc -if true then calculate loss fucntion
% The output is Dict, a struct with D,W,Z, Loss(the loss function value)

C=max(trlabels); % how many classes
[d, N]=size(X); % M is the data dimension, N is the # of samples
rng(opt.rng)

ind = randperm(N);
% D=X(:, ind(1:opt.K));
D = randn(d, opt.K);
Z=randn(opt.K,N);
W=randn(opt.K,opt.Q);
[M, ~, ~] = getMH1H2_t2(trlabels, Z);
Loss=zeros(1,opt.max_iter);


end % end of function file