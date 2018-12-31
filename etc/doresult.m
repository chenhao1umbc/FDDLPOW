result_beta = zeros(length(K),length(beta));
sparsity_beta = zeros(length(K),length(beta));
tr_sparsity_beta = zeros(length(K),length(beta));
result_betaWEEK = zeros(length(K),length(beta));

for ind1 = 1: length(beta)
    [opts]=loadoptions(K,lbmd,mu,Q,nu,beta(ind1), SNR);
    nm = ['SNR', num2str(SNR), opts.mixnm];
    if exist(nm, 'file')
        load(nm)                            
        % run prep_ZF 
        if exist('Dict')==1
            Dict_mix = Dict;
        end
        if cv == 1
            Z = sparsecoding_mix_cv(Dict_mix, Database, opts); %%%%% cv or test **************
        else
            Z = sparsecoding_mix_test(Dict_mix, Database, opts);
        end
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
        if cv == 1
            N_t = size(Database.cv_mixlabel, 2); %%%%% cv or test **************
        else 
            N_t = size(Database.test_mixlabel, 2); %%%%% cv or test **************
        end
        opts.ln_test = N_t/featln;
        opts.equal = pctrl.equal;
        [acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts);

        result_beta(ind1) = acc_all
        result_betaWEEK(ind1) = acc_weak_av
        sparsity_beta(ind1) = mean(sum(Z ~= 0))/K;
        tr_sparsity_beta(ind1) = mean(sum(Dict_mix.Z ~= 0))/K;

        [SW,SB]=calcfisher(Dict_mix.Z,Database.tr_label,opts);
        fWZ=trace(W'*SW*W)-trace(W'*SB*W)+norm(W'*Dict_mix.Z,'fro')^2;              

    end 
end