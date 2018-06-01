function Loss=DDLMD_Loss_mix_t2(X,trlabels,opt,W,D,Z)
% this fucntion is made to calculate the loss fucntion value
C=max(trlabels);
N = size(X,2);
Nc = N / C;
[SW,SB]=calcfisher(Z,trlabels,opt);
fWZ=trace(W'*SW*W)-trace(W'*SB*W)+norm(W'*Z,'fro')^2;

gWZ=0;
% indC=1:C;
% for ii=indC
%     ind=find(trlabels==ii);
%     mc=mean(Z(:,ind),2);
%     mc_not_ind=indC;
%     mc_not_ind(find(mc_not_ind==ii))=[];
%     for jj=mc_not_ind
%         ind=find(trlabels==jj);
%         mc_not=mean(Z(:,ind),2);
%         gWZ=gWZ+(mc'*W*W'*mc_not)^2;
%     end
% end
mc = Z*kron(eye(C),ones(Nc,1)/Nc);
WWt = W*W';
for c1=1:C
    for c2=1:C
        if c1 ~= c2
            gWZ = gWZ + (mc(:,c1)'*WWt*mc(:,c2))^2;
        end
    end
end

Loss=norm(X-D*Z,'fro')^2+opt.lambda1*sum(abs(Z(:)))+opt.mu*fWZ+opt.nu*gWZ;

end % end of the function file