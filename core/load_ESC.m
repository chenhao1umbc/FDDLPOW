function [database] = load_ESC(N_c, SNR, pctrl)

% This function is made to load testing, cross-validation, and traing data with their labels
% The dataset is called GRID, we only take 10 people out of 34. 
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
if nargin < 1
    N_c = 1;
    SNR = 2000;
    pctrl.equal = true; % equal power
elseif nargin < 2
    SNR = 2000;
    pctrl.equal = true; % equal power
elseif nargin < 3
    pctrl.equal = true; % equal power    
end

if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end

whichclass = 1: 10;% S1 to S10
nClass = length(whichclass);
trln = 20; % trainging 300 per class
cvln = 10; % cross-validation data is 50 per class
ttln = 10; % testing data is 40 per class
perClassln = cvln+ ttln + trln;
cvln_mix = 10; % cross-validation data is 100 mixture samples per combination
ttln_mix = 10; % testing data is 100 mixture samples per combination

%% loading data
% loading non-mixture data for training
% load 'sct_esc10_16_0.25_m2.mat'
% load 'sct_esc10_16_0.25_m3_renorm.mat'
load 'sct_esc10_16_0.25_m3_log.mat'
label = labels;% class index key


data = norm_data(data(:,:));
dwl = [data; label]; % data with labels

% rng(100) % run 20 different numbers to get the averaged result
ind = randperm(perClassln);
ind_tr = ind(1:trln);
ind_cv = ind(trln+1:trln+cvln);
ind_tt = ind(trln+cvln+1:end);

[d, ~] = size(dwl); % data dimension d 
n_tr_p = trln*featln;
n_cv_p = cvln*featln;
n_tt_p = ttln*featln;
n_tr = n_tr_p*nClass;
n_cv = n_cv_p*nClass;
n_tt = n_tt_p*nClass;

trls=zeros(1,n_tr);
cvls=zeros(1,n_cv);
dwl_tr=zeros(d,n_tr);
dwl_cv=zeros(d,n_cv);
dwl_tt=zeros(d,n_tt);

for ii = 1:nClass    
    dwl_tr(:,n_tr_p*(ii-1)+1:n_tr_p*ii) = dwl(:,(whichclass(ii)-1)*...
        featln*40+getindx(ind_tr,featln)); % training samples         
    dwl_cv(:,n_cv_p*(ii-1)+1:n_cv_p*ii) = dwl(:,(whichclass(ii)-1)*...
        featln*40+getindx(ind_cv,featln)); 
    dwl_tt(:,n_tt_p*(ii-1)+1:n_tt_p*ii) = dwl(:,(whichclass(ii)-1)*...
        featln*40+getindx(ind_tt,featln)); 
end

tr_dat = dwl_tr(1:d-1, :); % samples
cv_dat = dwl_cv(1:d-1, :);
tt_dat = dwl_tt(1:d-1, :);
trls = dwl_tr(end, :); % labels
cvls = dwl_cv(end, :);
ttls = dwl_tt(end, :);

if N_c  == 1
    cv_mixdat = cv_dat;
    cvmixls = cvls;
    tt_mixdat = tt_dat;
    ttmixls = ttls;
else 
    [cv_mixdat, cvmixls, tt_mixdat, ttmixls] = loadmixesc(pctrl);
end
database.pctrl = pctrl;
database.SNR = SNR;
database.N_c = N_c; % how many classes of signals mixed
database.featln = featln;
database.tr_data = tr_dat;
database.tr_label = trls;
database.cv_data = cv_dat;
database.cv_label = cvls;
database.test_data = tt_dat;
database.test_label = ttls;
database.cv_mixdata = cv_mixdat;
database.cv_mixlabel = cvmixls;
database.test_mixdata = tt_mixdat;
database.test_mixlabel = ttmixls;
database.cvln_mix  = cvln_mix; % samples for training in logistis regression

end % end of the function file