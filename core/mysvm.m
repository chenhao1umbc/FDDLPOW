function acc = mysvm(Xtr, Xcvortest, Database, cvortest, k)
% perform k-neareast neighbors

if nargin < 5
    k = 5;
end

if sum(cvortest) == 2
    error(' error from file myknn.m')
end

labels = Database.tr_label;
labels = aoos(Database.tr_label,Database.featln,size(Database.tr_label, 2));
if cvortest(1) % do cv
%     cvtlabel = Database.cv_label;
    cvtlabel = aoos(Database.cv_label,Database.featln,size(Database.cv_label, 2));
else
%     cvtlabel = Database.test_label;
    cvtlabel = aoos(Database.test_label,Database.featln,size(Database.test_label, 2));
end

Mdl = fitcecoc( Xtr',labels');
prelabel = predict(Mdl, Xcvortest');
a = prelabel - cvtlabel(1,:)';
acc = length(a(a ==0))/length(a);


end % end of the file


