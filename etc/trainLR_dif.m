function B = trainLR_dif(opt, Database)
% train logistic regression for mixture cases
% load data
SNR = opt.SNR;
featln = Database.featln;
cvln_mix = Database.cvln_mix;

tr_ln = size(Database.tr_data, 2);
cv_ln = size(Database.cv_data, 2);
trcv_ln = tr_ln + cv_ln;
P = size(Database.tr_data, 1);
N = tr_ln + cv_ln + featln*(15+20)*cvln_mix; % N_c <4
X_all = zeros(P, N); % P*N matrix
Y_all = zeros(N, 6); % 6 classes #################################

X_all(:, 1:(tr_ln + cv_ln)) = [Database.tr_data, Database.cv_data]; % samples
lb = [1,zeros(1,5)];
for counter =1:6
    Y_all(tr_ln/6*(counter-1)+1:tr_ln/6*counter, :) = ...
    kron(ones(tr_ln/6,1), circshift(lb, counter-1));
end
for counter =1:6
    Y_all(cv_ln/6*(counter-1)+1+tr_ln:cv_ln/6*counter+tr_ln, :) = ...
    kron(ones(cv_ln/6,1), circshift(lb, counter-1));
end
rt = './data/SNout_LMdata4/qn22/SNR_difpower/';
part1 = 'norm_mix';
part2 = ['db449_6classqn22_positive_renorm_snr',num2str(SNR),'power_'];
part3 = ['db449_6classqn22_negative_renorm_snr',num2str(SNR),'power_'];
counter = 7; %#####################################################

for N_c = 2:3 % N_c classes of mixtures
    Power = zeros(1,N_c); % only 0 0 0 
    c = combnk(1:6, N_c); % ble bt fhss1 zb
    for indCl=1:size(c,1)
        indClnm = c(indCl, :);
        nmdb1=[rt, part1,num2str(indClnm),part2,num2str(Power),'.mat']; 
        nmdb2=[rt, part1,num2str(indClnm),part3,num2str(Power),'.mat'];    
        load(nmdb1)
        load(nmdb2)
        % concatenate the positve and negative parts
        for ii=1:size(db.features,2)/featln
            db2.features(:,1+(ii-1)*featln:ii*featln)=...
                flip(flip(db2.features(:,1+(ii-1)*featln:ii*featln),2),1);
        end
        temp = [db.features; db2.features];
        X_all(:, trcv_ln+1+(counter-7)*cvln_mix*featln:trcv_ln+...
            (counter-6)*cvln_mix*featln) = temp(:, 1:featln*cvln_mix);
        lb = zeros(1, 6);
        lb(indClnm) = 1/length(indClnm);
        Y_all(trcv_ln+1+(counter-7)*cvln_mix*featln:trcv_ln+...
            (counter-6)*cvln_mix*featln, :) = kron(ones(50*featln,1), lb);
        counter = counter +1;
    end
end

load(opt.mixnm)
Database.test_mixdata = X_all;
Z_all = sparsecoding_mix_test(Dict_mix,Database,opts);

X = (Dict_mix.W'*Z_all)'; % projected 
Y = Y_all;

% train logistic regression
warning off
B = mnrfit(X, Y_all);
save(['SNR',num2str(SNR),'_beta', num2str(opt.beta), 'B_X_Y.mat'], 'B', 'X', 'Y');

end % end of the file