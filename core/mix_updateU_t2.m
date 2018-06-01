function U = mix_updateU_t2(W_in, M)
% Solves U update with constrain U^TU=I,
% input   trlabels is the training data labels
%         M is sparse coefficients mean per class, a matrix of K by C
%         W_in is the projection vector
% output is U, the orthogonal constrain

W = W_in;
Y = W'*M;        
[~,S,V] = svd(Y);
D = S'*S;
U = Y*V*diag(1./sqrt(diag(D)))*V';
    
end