function [D,Z,W,U,Loss,opt]=initdict_t2(X,trlabels,opt)
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

rng(opt.rng)

% check checking the existing Dictionary
% nm = ['esc_FDDLP',opt.dataset,'_k',num2str(opt.K),'_lmbd',num2str(opt.lambda1)...
%     ,'_mu',num2str(opt.mu),'_Q',num2str(opt.Q),'.mat' ];
fileexistance = exist(opt.Dictnm);
if fileexistance==2
    load(opt.Dictnm)
    D=Dict.D;
    Z=Dict.Z;
    W=Dict.W;
    M = getM_t2(opt.K, opt.C, opt.Nc, Z);
    U = mix_updateU_t2(W, M);        
    opt.max_iter= 120;% because of good initialization
    Loss=zeros(4,opt.max_iter);
    disp('good from  initdict_t2')
else    
    D=randn(opt.M_d,opt.K);
    Z=randn(opt.K,opt.N);
    W=randn(opt.K,opt.Q);
    M = getM_t2(opt.K, opt.C, opt.Nc, Z);
    U = mix_updateU_t2(W, M);
    Loss=zeros(4,opt.max_iter);
end 

end % end of function file