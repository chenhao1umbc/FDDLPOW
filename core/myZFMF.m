function acc =  myZFMF(labels_pre, Database, cvortest)
l = labels_pre(1,:);
if cvortest
    lt = Database.cv_label(1,1:Database.featln:end);
else
    lt = Database.test_label(1,1:Database.featln:end);
end
a = l - lt;
acc = sum(a ==0)/size(lt,2);



end