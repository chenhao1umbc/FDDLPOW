function [D, Z, W, U, V, Delta, Loss, opt]=initdict(X,trlabels,opt)
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
[M_d,N]=size(X); % M is the data dimension, N is the # of samples
rng(opt.rng)

% check checking the existing Dictionary
% nm=['DDLMD','_k',num2str(opt.K),'_lmbd',num2str(opt.lambda1)...
% ,'_mu',num2str(opt.mu),'_Q',num2str(opt.Q),'iter_',num2str(opt.max_iter),'.mat' ];
nm=['hello'];
fileexistance=exist(nm);
if fileexistance==2
    load(nm)
    D=Dict.D;
    Z=Dict.Z;
    W=Dict.W;
    [H1, ~, H3] = getMH1H2H3(trlabels, Z);
    Delta = eye(size(H3, 2));
    U = mix_updateU(W, Z, Delta, H3);
    V = mix_updateV(W, Z, H1);
    Loss=zeros(3,opt.max_iter);    
    opt.max_iter=50;% because of good initialization
else    
    D=randn(M_d,opt.K);
    Z=randn(opt.K,N);
    W=randn(opt.K,opt.Q);
    [H1, ~, H3] = getMH1H2H3(trlabels, Z);
    Delta = eye(size(H3, 2));
    U = mix_updateU(W, Z, Delta, H3);
    V = mix_updateV(W, Z, H1);
    Loss=zeros(3,opt.max_iter);
end 

end % end of function file