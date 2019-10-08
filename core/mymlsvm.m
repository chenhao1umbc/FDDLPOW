function [acc_weak, acc_weak_av, acc_all] = mymlsvm(Xtr, Xcvortest, cvortest, opts)
% perform k-neareast neighbors
if sum(cvortest) ~= 1
    error(' error from file myknn.m')
end

% reconstruct training labels
n_tr = size(Xtr, 2);
C = opts.C; % how many classes
labels = sum( kron(diag(1:C), ones(n_tr/C, 1)), 2);

% train svm model
Mdl = fitcecoc( Xtr',labels);
[~, scores] = predict(Mdl, Xcvortest');
[~, labels_pre] = sort(scores.', 1, 'descend');
[acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts);

end % end of the file
