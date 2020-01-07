function W = mix_updateW(opt, H_bar_i, S, H3, Delta, U, V,  Z)
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
beta = opt.beta;
Nc = opt.Nc;
M = Z*H3;

ZH_bar_iH_bar_iZt = 0;
ZHbar_iDeltaV = 0;
for ii = 1: opt.C
    ZH_bar_iH_bar_iZt = ZH_bar_iH_bar_iZt + Z(:, 1+ Nc*(ii-1): Nc*ii)*H_bar_i^2*Z(:, 1+ Nc*(ii-1): Nc*ii)';    
    ZHbar_iDeltaV = ZHbar_iDeltaV + Z(:, 1+ Nc*(ii-1): Nc*ii)*H_bar_i'*V{ii}'*Delta(ii);
end

W = (nu*M*M' + mu*Z*S*Z' +  beta* ZH_bar_iH_bar_iZt)\(nu*M*U'+ beta*ZHbar_iDeltaV);


end % end of fucntion file