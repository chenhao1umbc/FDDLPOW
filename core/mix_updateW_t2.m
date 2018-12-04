function W = mix_updateW_t2(opt, H1, H2, M, U, Z)
% this function is to update W with D and Z fixed, for the mixture cases
% input   X is the training data,a matrix M by N, N data samples
%         trlabels is the training data labels
%         D is Dictionary, a matrix of M by K, K atoms
%         Z is sparse coefficients, a matrix of K by N
%         opt is traing options
%             opt.K -the number of atoms in Dictionary
%             opt.Q -the projected dimensions
%             opt.lambda1 -lambda1 control the sparsity level
%             opt.mu -mu the coeffecient for fisher term
%             opt.nu - orthogonal term coefficient
%             opt.max_iter - the maximum iteration
%             opt.ploteig - plot eigenvalues if true
%             opt.losscalc -if true then calculate loss fucntion
% output is W, the projection matrix, K by Q

nu = opt.nu;
mu = opt.mu;
N = length(H1(:,1));
W = pinv(nu*M*M' + mu*Z*(eye(N)-H1)^2*Z' - mu*Z*(H1-H2)^2*Z' + mu*1.1*Z*Z')*nu*M*U';


end % end of fucntion file