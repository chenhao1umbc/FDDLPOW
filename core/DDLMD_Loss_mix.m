function Loss=DDLMD_Loss_mix(X,trlabels,opt,W,D,Z,U,V,delta)
% this fucntion is made to calculate the loss fucntion value

C = max(trlabels);
N = size(X,2);
Nc = N / C;
H1 = kron(eye(C),ones(Nc)/Nc);
H2 = ones(N)/N;
H3 = kron(eye(C),ones(Nc, 1)/Nc); % M = Z*H3
M1 = eye(N) - H1;
M2 = H1 - H2;
WtZ = W'*Z;
WtZM1 = WtZ*M1;
fWZ = norm(WtZM1, 'fro')^2 - norm(WtZ*M2, 'fro')^2 +norm(WtZ,'fro')^2;
gWZDelta = norm(WtZ*H3 -U, 'fro')^2;
cWZ = norm(WtZM1' - delta*V, 'fro')^2;

Loss=norm(X-D*Z,'fro')^2+opt.lambda1*sum(abs(Z(:)))+...
    opt.mu*fWZ+opt.nu*gWZDelta + opt.beta*cWZ;

end % end of the function file