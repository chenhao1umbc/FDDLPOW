function V = mix_updateV(W, Z, H1)
% Solves U update with constrain V^TV=I,
% input   H1 is squre matrix, kron(eye(C),ones(Nc)/Nc), C*N by C*N 
%         Z is sparse coefficients, K by N
%         W is the projection, K by Q
% output is V, the orthogonal constrain


[~, N]=size(Z);
M1 = eye(N) - H1;
Y = M1*Z'*W;

[u,~,v] = svd(Y);
V = u * eye(size(u,1), size(v, 1)) * v';    
end