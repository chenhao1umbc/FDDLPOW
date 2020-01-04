W = Dict.W;
C = max(Database.tr_label);
N = size(Database.tr_label,2);
Nc = N / C;
opts.C = C; % 6 classes
featln = Database.featln;
opts.n = Database.N_c;                
H0 = kron(eye(C),ones(Nc, 1)/Nc); % M = Z*H3
M = Dict.Z*H0;



if cvortest == 1
    N_t = size(Database.cv_mixlabel, 2); %%%%% cv or test **************
else 
    N_t = size(Database.test_mixlabel, 2); %%%%% cv or test **************
end
opts.ln_test = N_t/featln;
opts.equal = pctrl.equal;
opts.Ncombs = max(Database.cv_mixlabel);
