function acc = myknn(Xtr, Xcvortest, Database, cvortest, k)
% perform k-neareast neighbors

if nargin < 5
    k = 5;
end

if sum(cvortest) == 2
    error(' error from file myknn.m')
end

labels = Database.tr_label(1,:);
if cvortest(1) % do cv
    cvtlabel = Database.cv_label;
else
    cvtlabel = Database.test_label;
end

Mdl = fitcknn( Xtr',labels','NumNeighbors',k,'Standardize',1);
prelabel = predict(Mdl, Xcvortest');
a = prelabel - cvtlabel(1,:)';
acc = length(a(a ==0))/length(a);





end % end of the file
