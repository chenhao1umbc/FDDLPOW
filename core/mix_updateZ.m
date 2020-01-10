function Z = mix_updateZ(X,H_bar_i, H3, opt, W, D, Zin, U, V, Delta)
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
Nc = opt.Nc;
S = opt.S;
max_eig_S = opt.max_eig_S;
max_eig_H3 = opt.max_eig_H3;
H3H3t = opt.H3H3t;
H_bar_iSq = H_bar_i^2 ;
DtD = D'*D;
DtX = D'*X;
WWt = W*W';
WUH3t = W*U*H3';
normWWt = norm(WWt,'fro');

% calculate L, which is the Lipch. bound
L_term1 = norm(2*DtD,'fro');
L_term2 = 2 * opt.mu * normWWt * max_eig_S; 
L_term3 = 2*opt.nu*normWWt*max_eig_H3;  
L_term4 = 2* opt.beta * normWWt * sum_max_eig_Hhat_c;
L = L_term1 + L_term2 + L_term3 + L_term4;

% main loop
Z = fista(Z_in, L, opt.lambda1, opt, @calc_F, @grad_f);        

%% cost function
function cost = calc_F(Z_curr)   
    
    WtZ = W'*Z_curr;
    fisherterm=trace(WtZ*S*WtZ');

    % calc g(W, Z, U, Delta), orthogonal term
    M = Z_curr*H3;
    gwzu = norm(W'*M-U, 'fro')^2;            
    
    % calc Omega(W, Z), whitening term
    OmegaWZDeltaV = 0;    
    for ind = 1:opt.C
        OmegaWZDeltaV = OmegaWZDeltaV + ...
            norm(H_bar_i*WtZ(:, 1+Nc*(ind-1):Nc*ind)' - Delta(ind)*V{ind}, 'fro')^2;
    end
    
    % calc cost
    cost = norm(X - D*Z_curr,'fro')^2 + opt.lambda1*norm1(Z_curr) ...
        + opt.mu*fisherterm + opt.nu*gwzu + opt.beta*OmegaWZDeltaV;
end

%% gradiant of f
function g = grad_f(Z_curr)
    WWtZ = WWt*Z_curr;
    grad1 = DtD*Z_curr - DtX;
    grad2 = WWtZ*S;
    grad3 = WWtZ*H3H3t - WUH3t;
    
    grad4 = zeros(opt.K, opt.N);
    for ii = 1:opt.C    
        grad4(:, 1+Nc*(ii-1):Nc*ii)=  WWtZ(:, 1+Nc*(ii-1):Nc*ii)*H_bar_iSq - W*Delta(ii)*V{ii}*H_bar_i;
    end
    g = 2*grad1 + 2*opt.mu*grad2 + 2*opt.nu*grad3 + 2*opt.beta*grad4; 
end


end % end of this funciton file