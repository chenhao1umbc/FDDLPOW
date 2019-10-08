function Loss=DDLMD_Loss(X,trlabels,opt,W,D,Z)
% this fucntion is made to calculate the loss fucntion value

[SW,SB]=calcfisher(Z,trlabels,opt);
fWZ=trace(W'*SW*W)-trace(W'*SB*W)+ 1.1*norm(W'*Z,'fro')^2;
Loss=norm(X-D*Z,'fro')^2+opt.lambda1*sum(abs(Z(:)))+opt.mu*fWZ;

end % end of the function file