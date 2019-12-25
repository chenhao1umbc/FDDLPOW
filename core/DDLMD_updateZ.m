function Zout=DDLMD_updateZ(X,trlabels,opt,W,D,Zin)
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

C=max(trlabels);
N=length(trlabels);
Nc=N/C;
%diagONE=blockones(C,Nc);
H1 = kron(eye(C),ones(Nc)/Nc);
H2 = ones(N)/N;
M = (eye(N) - H1)^2 - (H1 - H2)^2 + 1.1*eye(N);
DtX = D'*X;
DtD = D'*D;
WWt = W*W';
opt.showprogress=opt.showconverge; % show the FISTA progress

% L = max(eig(2*D'*D)) + 4*opt.mu*max(eig(W*W'));
normWWt = norm(WWt,'fro');
L_term1 = 2*norm(DtD,'fro'); 
L_term2 = 2 * opt.mu * normWWt * norm(M,'fro');
L = L_term1 + L_term2 ;
sqrt_numelx = sqrt(numel(Zin));
for i = 1:opt.max_Ziter
    Zout=fista(Zin, L, opt.lambda1, opt, @calc_F, @grad);
    e = norm(Zin - Zout,'fro')/sqrt_numelx;
    if e < opt.threshold
        break;
    end
    Zin = Zout;
end
% convex function f
function cost = calc_F(Z)
    [SW,SB]=calcfisher(Z,trlabels,opt);
    fisherterm=trace(W'*SW*W)-trace(W'*SB*W)+ 1.1*norm(W'*Z, 'fro')^2;
    cost = norm((X - D*Z), 'fro')^2 + opt.mu*fisherterm + opt.lambda1*norm1(Z); 
end

% gradiant of f
function g = grad(Z)
    g = 2*(DtD*Z - DtX + opt.mu*WWt*Z*M);
end

end % end of the function file
