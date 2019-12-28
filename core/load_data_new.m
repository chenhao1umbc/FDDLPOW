function [database]=load_data_new(N_c, SNR, pctrl, rseed)

%This function is made to load the newly generated sacttering output with
%q[16, 0.05], n =2

% This function is made to load testing, cross-validation, and traing data with their labels
% output is database, struct with 
%         database.featln - scattering coefficient time lenght per sample
%         database.tr_data - training data 
%         database.tr_label - training labels
%         database.cv_data - cross validation data for non-mixture
%         database.cv_label - cross validation labels
%         database.cv_mixdata - cross-val data for mixture
%         database.cv_mixlabel - cross-val labels for mixture
%         database.test_mixdata - test data for mixture
%         database.test_mixlabel -test data labels for mixture
if nargin<1
    N_c = 2;
    SNR = 2000;
    pctrl.equal = 1; % equal power
elseif nargin<2
    SNR = 2000;
    pctrl.equal = 1; % equal power
elseif nargin<3
    pctrl.equal = 1; % equal power    
end
whichclass=[1:6];% bluetooth,BLE,wifi,zigbee
nClass=length(whichclass);
featln=8; % time length      <<<<*********************************** 
trln=300; % trainging 300 per class
cvln=74; % cross-validation data is 50 per class
ttln=75; % testing data is 40 per class
cvttln=cvln+ttln;
cvln_mix=100; % cross-validation data is 100 mixture samples per combination
ttln_mix=100; % testing data is 100 mixture samples per combination

%% loading data
% loading non-mixture data for training
nmdb1=['norm_db449_6classq16&0.05n2_positive_M2_snr', num2str(SNR) ,'.mat']; % positive
nmdb2=['norm_db449_6classq16&0.05n2_negative_M2_snr', num2str(SNR) ,'.mat']; % negative
load(nmdb1)
load(nmdb2)
% concatenate the positve and negative parts
for ii=1:449*6                
    db2.features(:,1+(ii-1)*featln:ii*featln)=...
        flip(flip(db2.features(:,1+(ii-1)*featln:ii*featln),2),1);
end
db.features=[db.features;db2.features];

rng(rseed);

ind = randperm(449);
ind_tr = ind(1:trln);
ind_cv = ind(trln+1:trln+cvln);
ind_tt = ind(trln+cvln+1:end);

n_r = size(db.features,1);
n_tr_p = trln*featln;
n_cv_p = cvln*featln;
n_tt_p = ttln*featln;
n_tr = n_tr_p*nClass;
n_cv = n_cv_p*nClass;
n_tt = n_tt_p*nClass;

trls=zeros(1,n_tr);
cvls=zeros(1,n_cv);
tr_dat=zeros(n_r,n_tr);
cv_dat=zeros(n_r,n_cv);
tt_dat=zeros(n_r,n_tt);

for ii=1:nClass  
    tr_dat(:,n_tr_p*(ii-1)+1:n_tr_p*ii)=db.features(:,(whichclass(ii)-1)*...
        featln*449+getindx(ind_tr,featln)); % training samples         
    cv_dat(:,n_cv_p*(ii-1)+1:n_cv_p*ii)=db.features(:,(whichclass(ii)-1)*...
        featln*449+getindx(ind_cv,featln)); 
    tt_dat(:,n_tt_p*(ii-1)+1:n_tt_p*ii)=db.features(:,(whichclass(ii)-1)*...
        featln*449+getindx(ind_tt,featln)); 
    
    trls(n_tr_p*(ii-1)+1:n_tr_p*ii)=ones(1,n_tr_p)*ii; %training labels
    cvls(n_cv_p*(ii-1)+1:n_cv_p*ii)=ones(1,n_cv_p)*ii; %cv labels
    ttls(n_tt_p*(ii-1)+1:n_tt_p*ii)=ones(1,n_tt_p)*ii;
end

if N_c ==1
    cv_mixdat = cv_dat;
    cvmixls = cvls;
    tt_mixdat = tt_dat;
    ttmixls = ttls;
else
    run ldd4_q16005n2
end
database.SNR = SNR;
database.N_c = N_c; % how many classes of signals mixed
database.featln=featln;
database.tr_data=tr_dat;
database.tr_label=trls;
database.cv_data=cv_dat;
database.cv_label=cvls;
database.test_data=tt_dat;
database.test_label=ttls;
database.cv_mixdata=cv_mixdat;
database.cv_mixlabel=cvmixls;
database.test_mixdata=tt_mixdat;
database.test_mixlabel=ttmixls;
database.cvln_mix =cvln_mix; % samples for training in logistis regression
end % end of the function file