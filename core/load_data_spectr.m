function [database]=load_data_spectr(N_c, SNR, pctrl)
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


trln=100; % trainging 300 per class
cvln=135; % cross-validation data is 50 per class
ttln=135; % testing data is 40 per class
cvttln=cvln+ttln;
trcvttln = cvttln+trln;
cvln_mix=50; % cross-validation data is 100 mixture samples per combination
ttln_mix=50; % testing data is 100 mixture samples per combination

%% data generation
load('/home/chenhao/Matlab/FDDLOW/data/total_dat');

rng(1)
ind = randperm(trcvttln);
ind_tr = ind(1:trln);
ind_cv = ind(trln+1:trln+cvln);
ind_tt = ind(trln+cvln+1:end);

featln = 199;
n_r = size(total_dat, 1);
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

whichclass=[1:nClass];
for ii=1:nClass  
    tr_dat(:,n_tr_p*(ii-1)+1:n_tr_p*ii)=total_dat(:,(whichclass(ii)-1)*...
        featln*trcvttln+getindx(ind_tr,featln)); % training samples         
    cv_dat(:,n_cv_p*(ii-1)+1:n_cv_p*ii)=total_dat(:,(whichclass(ii)-1)*...
        featln*trcvttln+getindx(ind_cv,featln)); 
    tt_dat(:,n_tt_p*(ii-1)+1:n_tt_p*ii)=total_dat(:,(whichclass(ii)-1)*...
        featln*trcvttln+getindx(ind_tt,featln)); 
    
    trls(n_tr_p*(ii-1)+1:n_tr_p*ii)=ones(1,n_tr_p)*ii; %training labels
    cvls(n_cv_p*(ii-1)+1:n_cv_p*ii)=ones(1,n_cv_p)*ii; %cv labels
    ttls(n_tt_p*(ii-1)+1:n_tt_p*ii)=ones(1,n_tt_p)*ii;
end

database.N_c = N_c; % how many classes of signals mixed
database.featln=featln;
database.tr_data=tr_dat;
database.tr_label=trls;
database.cv_data=cv_dat;
database.cv_label=cvls;
database.test_data=tt_dat;
database.test_label=ttls;
if N_c > 1
    cv_mixdat=[cv_mixdat,cv_dat_temp];
    cvmixls=[cvmixls,indCl*ones(1,size(cv_dat_temp,2))];
    database.cv_mixdata=cv_mixdat;
    database.cv_mixlabel=cvmixls;
    database.test_mixdata=tt_mixdat;
    database.test_mixlabel=ttmixls;
    database.cvln_mix =cvln_mix; % samples for training in logistis regression
end
end % end of the function file