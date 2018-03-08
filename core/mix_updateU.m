function U = mix_updateU(W, Z, Delta, H3)
% Solves U update with constrain U^TU=I,
% input   H3 is squre matrix, kron(eye(C),ones(Nc,1)/Nc), C*N by C
%         Z is sparse coefficients, K by N
%         Delta is C by C matrix
%         W is the projection, K by Q
% output is U, the orthogonal constrain

Y = W'*Z*H3*Delta';
      
[u0, ~, v0] = svd(Y);
U = u0 * eye(size(u0,1), size(v0, 1)) * v0';
    
end