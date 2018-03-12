function [acc, errortable, acc_av] = lr_test(Dict, Database, Z, B)

C = 6; % 6 classes
ln_test=length(Database.test_mixlabel);
true_labels = Database.test_mixlabel;
W = Dict.W;
wz = (W'*Z);
pre_prob = mnrval(B, wz');
[~,labels_pre] = sort(pre_prob, 2, 'descend');
labels_pre = labels_pre';

% run mix_n_ld34
run mix_ntable


end % end of the function file