function [acc_weak, acc_weak_av, acc_all] = mymlknn(Xtr, Xcvortest, cvortest, opts)
% perform k-neareast neighbors

if nargin < 5
    opts.k = 5;
    opts.C = 10;
end
if ~isfield(opts, 'k')
    opts.k = 5;
end
if sum(cvortest) ~= 1
    error(' error from file myknn.m')
end

labels_pre = mlknn(Xtr, Xcvortest,opts);
[acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts);

    function labels_pre = mlknn(train, cvtest, opts)
        % this function is calculated the nearest summation five samples for each class
        % return label vector for the test sample sorted by distance
        k = opts.k;
        C = opts.C;
        n_tr = size(train,2); % total training samples
        n_ts = size(cvtest,2); % total testing/cv samples
        dist = 1000*ones(n_tr, 1);
        result = ones(C, n_ts);
        ntrperC = n_tr/C;
        for ii = 1:n_ts
            for i = 1:n_tr            
                dist(i) = norm(train(:,i)-cvtest(:,ii));
            end
            for i = 1:C
                result(i,ii) = sum(mink(dist(1+(i-1)*ntrperC: ntrperC*i), k));
            end
        end    
        [~, labels_pre] = sort(result, 1, 'ascend');
    end % end of the function

end % end of the file
