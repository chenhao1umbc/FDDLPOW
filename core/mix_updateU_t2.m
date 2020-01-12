function U = mix_updateU_t2(W_in, M)
% Solves U update with constrain U^TU=I,
% input   trlabels is the training data labels
%         M is sparse coefficients mean per class, a matrix of K by C
%         W_in is the projection vector
% output is U, the orthogonal constrain

W = W_in;
Y = W'*M;        
[u,s,v] = svd(Y);
ds = diag(s);
if ds(end) ==0
    U = u*eye(size(s))*v';
else
    D = s'*s;
    U = Y*v*diag(1./sqrt(diag(D)))*v';
end
    
end