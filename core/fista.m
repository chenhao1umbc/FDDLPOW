function [X, iter] = fista(Xinit, L, lambda, opts, calc_F,grad)     
% * A Fast Iterative Shrinkage-Thresholding Algorithm for 
% Linear Inverse Problems.
% * Solve the problem: `X = arg min_X F(X) = f(X) + lambda||X||_1` where:
%   - `X`: variable, can be a matrix.
%   - `f(X)` is a smooth convex function with continuously differentiable 
%       with Lipschitz continuous gradient `L(f)` (Lipschitz constant of 
%       the gradient of `f`).
% * Syntax: `[X, iter] = FISTA(Xinit, L, lambda, opts, calc_F,grad)` where:
%   - INPUT:
%     + `grad`: a _function_ calculating gradient of `f(X)` given `X`.
%     + `Xinit`: initial guess.
%     + `L`: the Lipschitz constant of the gradient of `f(X)`.
%     + `lambda`: a regularization parameter, can be either a scalar or 
%       a weighted matrix.
%     + `opts`: a _structure_ variable describing the algorithm.
%       * `opts.max_iter`: maximum iterations of the algorithm. 
%           Default `300`.
%       * `opts.tol`: a tolerance, the algorithm will stop if difference 
%           between two successive `X` is smaller than this value. 
%           Default `1e-8`.
%       * `opts.show_progress`: showing `F(X)` after each iteration or not. 
%           Default `false`. 
%     + `calc_F`: optional, a _function_ calculating value of `F` at `X` 
%       via `feval(calc_F, X)`. 
% -------------------------------------
% Author: Tiep Vu, thv102, 4/6/2016
% (http://www.personal.psu.edu/thv102/)
% -------------------------------------
Linv = 1/L;    
lambdaLiv = lambda*Linv;
x_old = Xinit;
y_old = Xinit;
t_old = 1;
iter = 0;
cost=zeros(opts.max_iter,1);
E=zeros(opts.max_iter,1);
sqrt_numelx = sqrt(numel(Xinit));

%% MAIN LOOP
while  iter < opts.max_iter
    iter = iter + 1;
    x_new = shrinkage(y_old - Linv*feval(grad, y_old), lambdaLiv);
    t_new = 0.5*(1 + sqrt(1 + 4*t_old^2));
    y_new = x_new + (t_old - 1)/t_new * (x_new - x_old);
    
    % check stop criteria
    e = norm(x_new - x_old,'fro')/sqrt_numelx;
    if e < opts.threshold
        break;
    end
    
    % update
    x_old = x_new;
    t_old = t_new;
    y_old = y_new;
    
    % show progress
    if opts.showprogress        
        E(iter) = e;
        
        if opts.showcost
            cost(iter) = feval(calc_F, x_new);
        end
    end 
end
X = x_new;

%% show progress
if opts.showprogress
    figure(200);    

    if opts.showcost
    subplot(2,1,1)
    %hold on
    plot(cost(1:iter));
    legend('cost function value')
    title('cost function value updating Z, with W,D fixed')
    subplot(2,1,2)
    %hold on
    end
  
    semilogy(E(1:iter));
    legend('norm1(x_new - x_old)/numel(x_new)')
    title('Z difference, with W,D fixed')
    xlabel({'iterations';'-- from file "fista.m"'})
    pause(.1);
end

end
