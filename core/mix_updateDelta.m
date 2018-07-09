function Delta = mix_updateDelta(H_bar, Z, W, V, opt)
% Solves U update with constrain V^TV=I,
% input   M1 is H_bar
%         Z is sparse coefficients, K by N
%         W is the projection, K by Q
% output is Delta, 1 by C

Delta = ones(1, opt.C);
Nc = opt.Nc;
for ii = 1: opt.C
    H_bar_i = zeros(opt.N);
    H_bar_i(1+Nc*(ii-1):Nc*ii, 1+Nc*(ii-1):Nc*ii) =...
        H_bar(1+Nc*(ii-1):Nc*ii, 1+Nc*(ii-1):Nc*ii);% M1 is H_bar
    Delta(ii) = sum(sum(V{ii}.*(H_bar_i*Z*W)))/norm(V{ii}, 'fro')^2;
   
end



end % end of the file