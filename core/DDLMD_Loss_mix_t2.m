function Loss=DDLMD_Loss_mix_t2(X,trlabels,opt,W,D,Z, M, U)
% this fucntion is made to calculate the loss fucntion value
[SW,SB]=calcfisher(Z,trlabels,opt);
fWZ=trace(W'*SW*W)-trace(W'*SB*W)+1.1*norm(W'*Z,'fro')^2;

gWZ=0;
gWZ = norm(W'*M - U)^2;

Loss=norm(X-D*Z,'fro')^2+opt.lambda1*sum(abs(Z(:)))+opt.mu*fWZ+opt.nu*gWZ;

end % end of the function file