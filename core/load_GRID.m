function [database] = load_GRID(N_c, SNR, pctrl)

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
trln = 600; % trainging 300 per class
cvln = 200; % cross-validation data is 50 per class
ttln = 200; % testing data is 40 per class
perClassln = cvln+ ttln + trln;
cvln_mix = 100; % cross-validation data is 100 mixture samples per combination
ttln_mix = 100; % testing data is 100 mixture samples per combination

%% loading data
% loading non-mixture data for training
load 'GRID3k.mat' % data variable name is grid5k
label = zeros(2, 1e4); % 1000 per class, 10 classes
label(1, :) = sum(kron(diag(1:nClass), ones(1, perClassln)),1);% class index key
label(2, :) = 1:1e4; % primary key 0 to N for all the classes

grid3k = norm_data(grid3k); % normalize the data
dwl = [grid3k; label]; % data with labels

rng(100)
ind = randperm(perClassln);
ind_tr = ind(1:trln);
ind_cv = ind(trln+1:trln+cvln);
ind_tt = ind(trln+cvln+1:end);

[d, ~] = size(dwl); % data d
dwl_tr = zeros(d, trln*nClass);
dwl_cv = zeros(d, cvln*nClass);
dwl_tt = zeros(d, ttln*nClass);
for i = 1:nClass
    dwl_tr(:, trln*(i-1) + 1:trln*i ) = dwl(:, perClassln*(i-1) + ind_tr);
    dwl_cv(:, cvln*(i-1) + 1:cvln*i ) = dwl(:, perClassln*(i-1) + ind_cv);
    dwl_tt(:, ttln*(i-1) + 1:ttln*i ) = dwl(:, perClassln*(i-1) + ind_tt);
end
tr_dat = dwl_tr(1:d-2, :); % samples
cv_dat = dwl_cv(1:d-2, :);
tt_dat = dwl_tt(1:d-2, :);
trls = dwl_tr(d-1:end, :); % labels
cvls = dwl_cv(d-1:end, :);
ttls = dwl_tt(d-1:end, :);


if N_c  == 1
    cv_mixdat = cv_dat;
    cvmixls = cvls;
    tt_mixdat = tt_dat;
    ttmixls = ttls;
else % to be modified
    cv_mixdat = cv_dat;
    cvmixls = cvls;
    tt_mixdat = tt_dat;
    ttmixls = ttls;
end
database.SNR = SNR;
database.N_c = N_c; % how many classes of signals mixed
database.featln = 1;
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