function [D, Z, W, U, V, Delta, Loss, opt]=initdict(X, M1, H3, opt)
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

[M_d, ~]=size(X); % M is the data dimension, N is the # of samples
rng(opt.rng)

% check checking the existing Dictionary
nm = ['FDDLOW_mix','_k',num2str(opt.K),'_lmbd',num2str(opt.lambda1),...
    '_mu',num2str(opt.mu),'_Q',num2str(opt.Q),'_nu',num2str(opt.nu),...
    '_beta',num2str(-1),'.mat' ];
fileexistance=exist(nm);
if fileexistance==2
    load(nm)
    D=Dict_mix.D;
    Z=Dict_mix.Z;
    W=Dict_mix.W;
    Delta = ones(1, opt.C); 
    U = Dict_mix.U;   
    opt.max_iter=80;% because of good initialization
    Loss=zeros(3,opt.max_iter); 
else    
    D=randn(M_d,opt.K);
    Z=randn(opt.K,opt.N);
    W=randn(opt.K,opt.Q);d
    Delta = ones(1, opt.C); 
    U = mix_updateU(W, Z*H3);
    Loss=zeros(3,opt.max_iter);
end 
    V = mix_updateV(M1, Z, W, Delta, opt);
    
    
end % end of function file