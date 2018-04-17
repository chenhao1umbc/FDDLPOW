function [opts]=loadoptions(K,lambda1,mu,Q,nu,beta, SNR)
% this fucntion is made to load the options for dictionary learning
%             opt are the training/testing options with
%                 opt.K -the number of atoms in Dictionary
%                 opt.Q -the projected dimensions
%                 opt.lambda1 -lambda1 control the sparsity level
%                 opt.mu -mu the coeffecient for fisher term
%                 opt.max_iter - the maximum iteration
%                 opt.losscalc -if true then calculate loss fucntion
%                 opt.k_neighors - how many neighors
%                 opt.kNNmethod -1 for averaging; 2 for max pooling needs more computation
if nargin < 5
    nu=1000;
    beta = 1;
    SNR = 2000;
end
opts.SNR = SNR;
opts.K=K;
opts.Q=Q;
opts.lambda1=lambda1;                                        
opts.mu=mu;
opts.showconverge = false;
%opts.threshold=1e-8;% for FISTA sparse coding
opts.savedict=true;
opts.max_iter=150;
opts.losscalc=true;
opts.k_neighbors=3; % how many neighors
opts.kNNmethod=1; % 1 for averaging; 2 for max pooling needs more computation
opts.dataset = '4'; % LMdata 4
opts.Dictnm =['FDDLOW','_k',num2str(opts.K),'_lmbd',num2str(opts.lambda1)...
    ,'_mu',num2str(opts.mu),'_Q',num2str(opts.Q),'.mat' ];

% for mixture algorithm
opts.nu=nu;
opts.beta = beta;
opts.mixcase=true;
opts.mixnm =['FDDLOW_mix','_k',num2str(opts.K),'_lmbd',num2str(opts.lambda1)...
    ,'_mu',num2str(opts.mu),'_Q',num2str(opts.Q),'_nu',num2str(opts.nu),...
    '_beta',num2str(opts.beta),'.mat' ];

end % end of the function file