function [acc_weak, acc_weak_av, acc_all] = mymlknn(Xtr, Xcvortest, cvortest, opts)
% perform k-neareast neighbors

if nargin < 5
    opts.k = 3;
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

end % end of the file
