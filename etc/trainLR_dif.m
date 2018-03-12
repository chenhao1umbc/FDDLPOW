% train logistic regression for mixture cases

%% load data
SNR = 2000;
featln = 4;

P = size(Database.tr_data, 1);
N = size(Database.tr_data, 2)+size(Database.cv_data, 2)+featln*
449*4*6+(15+20)*4*50
X_all = zeros(P, N); % P*N matrix
Y_all = zeros(N, 6); % 6 classes

nmdb1=['norm_db449_6classqn22_positive_renorm_snr', num2str(SNR),'.mat']; % positive
nmdb2=['norm_db449_6classqn22_negative_renorm_snr', num2str(SNR),'.mat']; % negative
load(nmdb1)
load(nmdb2)

% concatenate the positve and negative parts
for ii=1:449*6                
    db2.features(:,1+(ii-1)*featln:ii*featln)=...
        flip(flip(db2.features(:,1+(ii-1)*featln:ii*featln),2),1);
end
X_all(:, 1:449*6*featln) = [db.features;db2.features]; % samples
lb = [1,zeros(1,5)];
for counter =1:6
    Y_all(featln*449*(counter-1)+1:featln*449*counter, :) = ...
    kron(ones(449*featln,1), circshift(lb, counter-1));
end
rt = './data/SNout_LMdata4/qn22/SNR_difpower/';
part1 = 'norm_mix';
part2 = ['db449_6classqn22_positive_renorm_snr',num2str(SNR),'power_'];
part3 = ['db449_6classqn22_negative_renorm_snr',num2str(SNR),'power_'];
counter = 7;

for N_c = 2:3
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
        X_all(:, 10777+(counter-7)*50*featln:10776+(counter-6)*50*featln) = temp(:, 1:featln*50);
        lb = zeros(1, 6);
        lb(indClnm) = 1/length(indClnm);
        Y_all(10777+(counter-7)*50*featln:10776+(counter-6)*50*featln, :) = ...
            kron(ones(50*featln,1), lb);
        counter = counter +1;
    end
end

load(['SNR',num2str(SNR),'DDLMDmix4_k100_lmbd1.5_mu0.1_Q16_nu1000iter_100.mat'])
Database.test_mixdata = X_all;
Z_all = sparsecoding_mix_test(Dict_mix,Database,opts);

X = (Dict_mix.W'*Z_all)'; % projected 
Y = Y_all;

% train logistic regression
warning off
B = mnrfit(X, Y_all);
save(['SNR',num2str(SNR),'B_X_Y_mix.mat'], 'B', 'X', 'Y');


toc