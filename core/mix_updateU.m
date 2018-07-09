function U = mix_updateU(W, Z, H3)
% Solves U update with constrain U^TU=I,
% input   H3 is squre matrix, kron(eye(C),ones(Nc,1)/Nc), C*N by C
%         Z is sparse coefficients, K by N
%         Delta is C by C matrix
%         W is the projection, K by Q
% output is U, the orthogonal constrain

Y = W'*Z*H3;        
[~,S,V] = svd(Y);
D = S'*S;
U = Y*V*diag(1./sqrt(diag(D)))*V';

    
end