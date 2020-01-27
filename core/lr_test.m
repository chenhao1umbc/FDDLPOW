function [acc_weak, acc_weak_av, acc_all] = lr_test(Dict, Database, Z, B, pctrl)

opts.C = max(Database.tr_label); % 6 classes
featln = Database.featln;
opts.n = Database.N_c;
opts.Ncombs = max(Database.cv_mixlabel);
N = size(Database.test_mixlabel, 2);
opts.ln_test = N/featln;
opts.equal = pctrl.equal;

W = Dict.W;
wz = W'*Z; %(W'*aoos(Z, featln, N));
pre_prob = mnrval(B, wz');
[~,labels_pre] = sort(pre_prob, 2, 'descend');
labels_pre = labels_pre';

if pctrl.if2weak == 0
    [acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts);
else
    [acc_weak, acc_weak_av, acc_all] = calc_labels2w(labels_pre, opts);
end

end % end of the function file
