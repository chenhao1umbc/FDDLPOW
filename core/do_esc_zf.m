W = Dict.W;
C = max(Database.tr_label);
N = size(Database.tr_label,2);
Nc = N / C;
opts.C = C; % 10 classes
featln = Database.featln;
opts.n = Database.N_c;                
H3 = kron(eye(C),ones(Nc, 1)/Nc); % M = Z*H3
M = Dict.Z*H3;
% zero forcing
H = W'*M;
result = pinv(H)*W'*Z;
[~, labels_pre] = sort(result, 1, 'descend');
opts.Ncombs = max(Database.cv_mixlabel);
opts.ln_test = size(Database.test_mixlabel, 2)/featln;
opts.equal = pctrl.equal;
[acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts);