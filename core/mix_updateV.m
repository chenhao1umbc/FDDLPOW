function [V, WtZHhatc] = mix_updateV(H_bar_i, Z, W, Delta, opt)
% Solves U update with constrain V^TV=I,
% input   H1 is squre matrix, kron(eye(C),ones(Nc)/Nc), C*N by C*N 
%         Z is sparse coefficients, K by N
%         W is the projection, K by Q
%         Delta is 1 by C
% output is V, the orthogonal constrain, is a cell 1 by C, each cell is
%           Nc by Nc

Nc = opt.Nc;
V = cell(1, opt.C);
WtZHhatc = cell(1, opt.C);
for ii = 1: opt.C
    Y = H_bar_i*Z(:, 1+ Nc*(ii-1): Nc*ii)'*W/Delta(ii);
    WtZHhatc{ii} = Y'*Delta(ii);
    [u,~,v] = svd(Y);
    V{ii} = v * eye(size(v,2), size(u, 2)) * u';    
end

end % end of the file