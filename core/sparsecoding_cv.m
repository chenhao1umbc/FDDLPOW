function Z_cv=sparsecoding_cv(Dict,database,opt)
% this function is used for sparse coding for the testing/cv data
% min norm2(x-Dz)^2+lambda1*norm1(z)

D=Dict.D;
X=database.cv_data;
opt.max_iter = 500;
opt.threshold=1e-4;
opt.showprogress=false; % show the FISTA progress
rng(0)
Zinit=randn(opt.K,size(X,2));

L = max(eig(2*D'*D));  
Z_cv=fista(Zinit, L, opt.lambda1, opt, @calc_F,@grad);

% cost function
function cost = calc_F(Z)
    cost = calc_f(Z) + opt.lambda1*norm1(Z);
end 
% convex function f
function cost = calc_f(Z)        
    cost =norm((X - D*Z), 'fro')^2; 
end 
% gradiant of f
function g = grad(Z)
    g= 2*D'*(D*Z-X);
end

end % end of the function file