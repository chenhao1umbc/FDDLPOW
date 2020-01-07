function Loss=Loss_mix(X, H_bar_i, H3, S, opt,W,D,Z,U,V,Delta)
% this fucntion is made to calculate the loss fucntion value

WtZ = W'*Z;
% WtZM1 = WtZ*M1;
% fWZ = norm(WtZM1, 'fro')^2 - norm(WtZ*M2, 'fro')^2 + 1.1*norm(WtZ,'fro')^2;
fWZ = trace(WtZ*S*WtZ');
gWZDelta = norm(WtZ*H3 -U, 'fro')^2;

OmegaWZDeltaV = 0;
Nc = opt.Nc;
for ii = 1:opt.C    
    OmegaWZDeltaV = OmegaWZDeltaV + ...
        norm(WtZ(:, 1+ Nc*(ii-1): Nc*ii)*H_bar_i - Delta(ii)*V{ii}, 'fro')^2;
end

Loss=norm(X-D*Z,'fro')^2+opt.lambda1*sum(abs(Z(:)))+...
    opt.mu*fWZ+opt.nu*gWZDelta + opt.beta*OmegaWZDeltaV;

end % end of the function file