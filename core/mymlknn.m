function [acc_weak, acc_weak_av, acc_all] = mymlknn(Xtr, Xcvortest, Database, cvortest, k)
% perform k-neareast neighbors

if nargin < 5
    k = 5;
end
if sum(cvortest) ~= 1
    error(' error from file myknn.m')
end

opts.C = C; % 10 classes
featln = Database.featln;
opts.n = Database.N_c; 
opts.Ncombs = max(Database.cv_mixlabel);
opts.equal = pctrl.equal;

if cvortest(1) % do cv
    opts.ln_test = size(Database.cv_mixlabel, 2)/featln;
else
    opts.ln_test = size(Database.test_mixlabel, 2)/featln;
end

labels_pre = mlknn(Xtr, Xcvortest, k);
[acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts);

    function labels_pre = mlknn(train, cvtest, k)
        % this function is calculated the nearest summation five samples for each class
        % return label vector for the test sample sorted by distance
        n_tr = size(train,2); % total training samples
        n_ts = size(cvtest,2); % total testing/cv samples
        dist = 1000*ones(n_tr, 1);
        for ii = 1:n_ts
            for i = 1:n_tr            
                dist(i) = norm(train(:,i),cvtest(:,ii));
            end
            for i = 1:C
                result(i, = mink(dist((i-1)*), k);
            end
        end    
        labels_pre
    end % end of the function

end % end of the file
