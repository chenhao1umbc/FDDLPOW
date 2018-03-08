function D=DDLMD_updateD(X,opt,D,Z)
% this function is to update W with D and Z fixed
% input   X is the training data,a matrix M by N, N data samples
%         D is Dictionary, a matrix of M by K, K atoms
%         Z is sparse coefficients, a matrix of K by N
%         opt is traing options
%             opt.K -the number of atoms in Dictionary
%             opt.Q -the projected dimensions
%             opt.lambda1 -lambda1 control the sparsity level
%             opt.mu -mu the coeffecient for fisher term
%             opt.max_iter - the maximum iteration
%             opt.threshold - threshold to stop the iteration
%             opt.showconverge - show norm(D_old-D_new) if true
%             opt.losscalc -if true then calculate loss fucntion
% output is D, the dictionary, matrix M by K
% the algorithm is from J. Mairal, F. Bach, J. Ponce and Guillermo Sapiro,"Online Learning for Matrix 
% Factorization and Sparse"

[K,N]=size(Z);
%A=0;
%B=0;
% calc A and B 
%for ii=1:N
%    A=A+Z(:,ii)*Z(:,ii)'; % A should be K by K matrix
%    B=B+X(:,ii)*Z(:,ii)';% B should be M by K matrix
%end
A = Z*Z';
B = X*Z';

% main iteration to get D
dist = zeros(1,opt.max_iter);
for ii=1:opt.max_iter
    D_old=D;
    % loop through each atom
    for jj=1:K
%         % tackle Ajj is 0
%         if A(jj,jj)==0
%             A(jj,jj)=A(jj,jj)+1e-10;
%         end
        uj=(B(:,jj)-D*A(:,jj))/(A(jj,jj)+1e-10)+D(:,jj);
        D(:,jj)=uj/max([norm(uj),1]);
    end        
    dist(ii)=norm(D-D_old,'fro');
    if dist(ii)<opt.threshold
        break
    end
end

% plot convergence
if opt.showconverge
    figure(100);
    %hold on
    semilogy(dist);
    title('norm(D\_old-D\_new)')
    xlabel({'iterations';'-- from file "DDLMD\_updateD.m"'})
    pause(.1);
end

end % end of the function file

