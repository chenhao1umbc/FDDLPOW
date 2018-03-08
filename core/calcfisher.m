function [SW,SB]=calcfisher(Z,trlabels,opt)

SW=zeros(opt.K);
SB=zeros(opt.K);
C=max(trlabels); % number of classes
m=mean(Z,2);
% sum over each class
for ii=1:C
    ind=find(trlabels==ii);
    mc=mean(Z(:,ind),2);
    Nc=length(ind);
    Zc=Z(:,ind);
    % sum over within class
    for jj=1:Nc        
        SW=SW+(Zc(:,jj)-mc)*(Zc(:,jj)-mc)';
    end
    SB=SB+Nc*(m-mc)*(m-mc)';
end

end % end of function file