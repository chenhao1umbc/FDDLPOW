function acc = myknn(Xtr, Xcvortest, Database, cvortest, k)
% perform k-neareast neighbors

if nargin < 5
    k = 5;
end

if sum(cvortest) ~= 1
    error(' error from file myknn.m')
end

labels = Database.tr_label;
% labels = aoos(Database.tr_label,Database.featln,size(Database.tr_label, 2));
if cvortest(1) % do cv
%     cvtlabel = Database.cv_label; % only used when DSS+knn
    cvtlabel = aoos(Database.cv_label,Database.featln,size(Database.cv_label, 2));
else
%     cvtlabel = Database.test_label; % only used when DSS+knn
    cvtlabel = aoos(Database.test_label,Database.featln,size(Database.test_label, 2));
end

Mdl = fitcknn( Xtr',labels','NumNeighbors',k,'Standardize',1);
prelabel = predict(Mdl, Xcvortest');

a = prelabel - cvtlabel(1,:)';
acc = length(a(a ==0))/length(a);





end % end of the file
