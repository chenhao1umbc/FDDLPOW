function V = mix_updateV(Y_in, delta)
% Solves U update with constrain V^TV=I,
% input   H1 is squre matrix, kron(eye(C),ones(Nc)/Nc), C*N by C*N 
%         Z is sparse coefficients, K by N
%         W is the projection, K by Q
% output is V, the orthogonal constrain

Y = Y_in/delta;
[u,~,v] = svd(Y);
V = u * eye(size(u,1), size(v, 1)) * v';    

end