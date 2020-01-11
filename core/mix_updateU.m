function U = mix_updateU(W, Z, H3)
% Solves U update with constrain U^TU=I,
% input   H3 is squre matrix, kron(eye(C),ones(Nc,1)/Nc), C*N by C
%         Z is sparse coefficients, K by N
%         Delta is C by C matrix
%         W is the projection, K by Q
% output is U, the orthogonal constrain

Y = W'*Z*H3;        
[u,s,v] = svd(Y);
if s(end) ==0
    U = u*eye(size(s))*v';
else
    D = s'*s;
    U = Y*v*diag(1./sqrt(diag(D)))*v';
end
    
end