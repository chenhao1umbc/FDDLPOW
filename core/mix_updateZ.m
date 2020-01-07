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
M2M2t = opt.M2M2t;
M1M1t = opt.M1M1t;
H3H3t = opt.H3H3t;
H_bar_iSq = H_bar_i^2 ;
DtD = D'*D;
DtX = D'*X;
WWt = W*W';
WUH3t = W*U*H3';
normWWt = norm(WWt,'fro');

% calculate L, which is the Lipch. bound
L_term1 = norm(2*DtD,'fro');
L_term2 = 2 * opt.mu * normWWt * norm(M1M1t-M2M2t+eye(opt.N),'fro'); 
L_term3 = 2*opt.nu*normWWt*norm(H3H3t, 'fro');  

sum_norm_H_bar_i = 0;
for i = 1: opt.C
    sum_norm_H_bar_i = sum_norm_H_bar_i + norm(H_bar_iSq, 'fro');
end
L_term4 = 2* opt.beta * normWWt * sum_norm_H_bar_i;
L = L_term1 + L_term2 + L_term3 + L_term4;

% main loop
dZ = inf(opt.max_Ziter, 1);
cost = dZ;
for iter = 1:opt.max_Ziter     
    Z_old = Z;       
    Z = fista(Z_old, L, opt.lambda1, opt, @calc_F, @grad_f);        
    
    % check convergence
    dZ(iter) = norm(Z - Z_old,'fro')/sqrt(numel(Z));    
    if dZ(iter) < opt.Zthreshold
        break;
    end
    
    if opt.showconverge
        iter
        cost(iter) = calc_F(Z_old);
        figure(210);
        subplot(2,1,1);
        plot(cost(1:iter));
        title('Cost function');
        
        subplot(2,1,2);
        plot(dZ(1:iter));
        title('||Z - Z_old||_F/sqrt(KN)');
        xlabel({'Iterations';'--from mix\_updateZ.m'})
        pause(.1);
    end
end

%% cost function
function cost = calc_F(Z_curr)   
    
    WtZ = W'*Z_curr;
    % calc fisherterm
    SW = Z_curr*M1M1t*Z_curr';
    SB = Z_curr*M2M2t'*Z_curr';
    fisherterm=trace(W'*SW*W)-trace(W'*SB*W)+1.1*norm(W'*Z_curr, 'fro')^2;

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
    WWtZM1M1t = WWtZ*M1M1t;
    grad1 = DtD*Z_curr - DtX;
    grad2 = WWtZM1M1t - WWtZ*M2M2t + WWtZ;
    grad3 = WWtZ*H3H3t - WUH3t;
    
    grad4 = zeros(opt.K, opt.N);
    for ii = 1:opt.C    
        grad4(:, 1+Nc*(ii-1):Nc*ii)=  WWtZ(:, 1+Nc*(ii-1):Nc*ii)*H_bar_iSq - W*Delta(ii)*V{ii}*H_bar_i;
    end
    g = 2*grad1 + 2*opt.mu*grad2 + 2*opt.nu*grad3 + 2*opt.beta*grad4; 
end


end % end of this funciton file