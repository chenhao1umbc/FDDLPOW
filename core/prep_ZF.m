
if exist('Dict')==1
    Dict_mix = Dict;
end
Z = sparsecoding_mix_test(Dict_mix, Database, opts);
W = Dict_mix.W;
C = max(Database.tr_label);
N = size(Database.tr_label,2);
Nc = N / C;
opts.C = C; % 6 classes
featln = Database.featln;
opts.n = Database.N_c;
H3 = kron(eye(C),ones(Nc, 1)/Nc); % M = Z*H3
M = Dict_mix.Z*H3;
% zero forcing
H = W'*M;
result = pinv(H)*W'*aoos(Z,featln,size(Z, 2));
[~, labels_pre] = sort(result, 1, 'descend');


opts.Ncombs = max(Database.cv_mixlabel);
N_t = size(Database.test_mixlabel, 2); % test signal length
opts.ln_test = N_t/featln;
opts.equal = pctrl.equal;