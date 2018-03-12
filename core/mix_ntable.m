% to see per combination accuracy
n = Database.N_c;
Ncombs=max(Database.cv_mixlabel);
c = combnk(1:6,n); % ble bt fhss1 zb
acc=zeros(C,1); % only score for true-positive and true negative
acc0=zeros(C,ln_test/Ncombs); % only score for true-positive and true negative
for indCl=1:Ncombs
    indClnm(1:n) = c(indCl,:);
    indClnm(n+1:6) = Theother3(indClnm(1:n));
    for ii=1:ln_test/Ncombs
        for jj = 1:n
            acc0(indClnm(jj),ii)=length(find(labels_pre(1:n,ii+...
                (indCl-1)*ln_test/Ncombs)==indClnm(jj)));
            acc(indClnm(jj))=acc(indClnm(jj))+length(find...
                (labels_pre(1:n,ii+(indCl-1)*ln_test/Ncombs)==indClnm(jj)));
        end
        for jj = n+1:6
            acc0(indClnm(jj),ii)=length(find(labels_pre(n+1:6,ii+...
                (indCl-1)*ln_test/Ncombs)==indClnm(jj)));
            acc(indClnm(jj))=acc(indClnm(jj))+length(find...
                (labels_pre(n+1:6,ii+(indCl-1)*ln_test/Ncombs)==indClnm(jj)));
        end
    end 
    errortable(:,indCl)=mean(acc0,2);
end
acc=acc/ln_test;

errortable = 1 - errortable;

% overall accuracy
if n==1
    acc_av = 1 - sum(diag(errortable))/6;
else
    acc_av=sum(acc)/C;
end