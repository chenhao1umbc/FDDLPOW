function W=DDLMD_updateW(trlabels,opt,Z)
% this function is to update W with D and Z fixed
% input   X is the training data,a matrix M by N, N data samples
%         trlabels is the training data labels
%         D is Dictionary, a matrix of M by K, K atoms
%         Z is sparse coefficients, a matrix of K by N
%         opt is traing options
%             opt.K -the number of atoms in Dictionary
%             opt.Q -the projected dimensions
%             opt.lambda1 -lambda1 control the sparsity level
%             opt.mu -mu the coeffecient for fisher term
%             opt.max_iter - the maximum iteration
%             opt.ploteig - plot eigenvalues if true
%             opt.losscalc -if true then calculate loss fucntion
% output is W, the projection matrix, K by Q

[SW,SB]=calcfisher(Z,trlabels,opt);
A = SW-SB+1.1*Z*Z';
A = (A +A')/2;
[V,D]=eig(A);  % lbd is eig value, increasing to right
[d_sorted,I] = sort(diag(D));
ind = find(d_sorted> 1e-7); % in case A is low rank
Wdimension = min(opt.Q, length(ind));% which eigen vector to pick
W = V(:,I(ind(1:Wdimension)));

if opt.ploteig
    figure(300);
    %hold on
    semilogy(d_sorted);
    title('eignvalues of SW-SB+Z*Z^T')
    xlabel('--from file DDLMD\_updateW.m')
    pause(.1);
end

end % end of fucntion file