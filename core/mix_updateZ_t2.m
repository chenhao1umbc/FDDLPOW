function Z = mix_updateZ_t2(X, trlabels, opt, W, D, Zin, U,H3, S, H3H3t,max_eig_S, max_eig_H3)
% this function is to update W with D and Z fixed
% input   X is the training data,a matrix M by N, N data samples
%         trlabels is the training data labels
%         W is projection, a matrix of K by N
%         D is Dictionary, a matrix of M by K, K atoms
%         opt is traing options
%             opt.K -the number of atoms in Dictionary
%             opt.Q -the projected dimensions
%             opt.lambda1 -lambda1 control the sparsity level
%             opt.mu -mu the coeffecient for fisher term
%             opt.max_iter - the maximum iteration
%             opt.threshold - threshold to stop the iteration
%             opt.losscalc -if true then calculate loss fucntion
% output is Z, the sparse coefficient,a matrix K by N
% the FISTA algorithm is from A. Beck and M. Teboulle, "A fast iterative shrinkage-thresholding
% algorithm for linear inverse problems", SIAM Journal on Imaging Sciences, vol. 2, no. 1, pp. 183�202, 2009

% pre-cacularion

Z = Zin;
DtD = D'*D;
DtX = D'*X;
WWt = W*W';
WUH3t = W*U*H3';
% normWWt = norm(WWt,'fro');
% L_term1 = norm(2*DtD,'fro'); 
% L_term2 = 2 * opt.mu * normWWt * norm(M1M1t-M2M2t+1.1*eye(N),'fro');
% L_term3 = 2*opt.nu*normWWt*norm(H3H3t, 'fro'); 
eigWWt = max(eig(WWt));
L_term1 = max(eig(2*DtD)); 
L_term2 = 2 * opt.mu * eigWWt * max_eig_S;
L_term3 = 2*opt.nu*eigWWt*max_eig_H3; 
L = L_term1 + L_term2 + L_term3;

% main loop       
Z = fista(Z, L, opt.lambda1, opt, @calc_F, @grad_f);        


%% cost function
function cost = calc_F(Z_curr)
    % calc fisherterm

    [SW,SB]=calcfisher(Z_curr,trlabels,opt);
    fisherterm=trace(W'*SW*W)-trace(W'*SB*W)+1.1*norm(W'*Z_curr, 'fro')^2;

    % calc g(W, Z, U)
    M = Z_curr*H3;
    gwzu = norm(W'*M-U, 'fro');            
    
    % calc cost
    cost = norm(X - D*Z_curr,'fro')^2 + opt.lambda1*norm1(Z_curr) ...
        + opt.mu*fisherterm + opt.nu*gwzu;
end

%% gradiant of f
function g = grad_f(Z_curr)
    grad1 = DtD*Z_curr - DtX;
    grad2 = WWt* Z_curr*S;
    grad3 = WWt*Z_curr*H3H3t - WUH3t;
    g = 2*grad1 + 2*opt.mu*grad2 + 2*opt.nu*grad3; 
end


end % end of this funciton file