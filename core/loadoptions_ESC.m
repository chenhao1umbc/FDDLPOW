function [opts]=loadoptions_ESC(whichtable, K,lambda1,mu,Q,nu,beta, SNR)
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
if nargin < 6
    nu=1000;
    beta = 1;
    SNR = 2000;
end
opts.C = 10; % total of C = 6 classes
opts.rng = 0; % for dictionary initialization
opts.SNR = SNR;
opts.K=K;
opts.Q=Q;
opts.lambda1=lambda1;                                        
opts.mu=mu;
opts.showconverge = false;
%opts.threshold=1e-8;% for FISTA sparse coding
opts.savedict=true;
opts.max_iter=100;
opts.losscalc=true;
opts.k_neighbors=5; % how many neighors
opts.kNNmethod=1; % 1 for averaging; 2 for max pooling needs more computation
opts.dataset = 'ESC10'; % LMdata 4

% for mixture algorithm
opts.nu=nu;
opts.beta = beta;
opts.mixcase=true;

if whichtable == 1
opts.Dictnm =['esc_FDDLP',opts.dataset,'_k',num2str(opts.K),'_lmbd',num2str(opts.lambda1)...
    ,'_mu',num2str(opts.mu),'_Q',num2str(opts.Q),'.mat' ];
end

if whichtable == 2
opts.Dictnm =['esc_FDDLPO',opts.dataset,'_k',num2str(opts.K),'_lmbd',num2str(opts.lambda1)...
    ,'_mu',num2str(opts.mu),'_nu',num2str(opts.nu), '_Q',num2str(opts.Q),'.mat' ];
end

if whichtable == 3
opts.Dictnm =['esc_FDDLPOW',opts.dataset,'_k',num2str(opts.K),'_lmbd',num2str(opts.lambda1),...
    '_mu',num2str(opts.mu),'_nu',num2str(opts.nu),...
    '_beta',num2str(opts.beta), '_Q',num2str(opts.Q),'.mat' ];
end


end % end of the function file