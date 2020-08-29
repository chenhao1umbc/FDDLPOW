function Delta = mix_updateDelta(WtZHhatc, V, opt)
% Solves U update with constrain V^TV=I,
% input   M1 is H_bar
%         Z is sparse coefficients, K by N
%         W is the projection, K by Q
% output is Delta, 1 by C

Delta = ones(1, opt.C);
for ii = 1: opt.C
%     Delta(ii) = sum(sum(V{ii}.*(H_bar_i*Z(:, 1+ Nc*(ii-1): Nc*ii)'*W)))/norm(V{ii}, 'fro')^2;   
Delta(ii) = trace(V{ii}'*WtZHhatc{ii})/norm(V{ii}, 'fro')^2;
end

end % end of the file