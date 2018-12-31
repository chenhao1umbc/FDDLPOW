function [cv_mixdat, cvmixls, tt_mixdat, ttmixls] = loadmixesc(pctrl, N_c, SNR )

% input N_c is how many classes in the mixture 
% SNR is the signal to noise ratio
% pctrl is the power of each components in the mixture, 0-0 (equal power),
% 3db 0-3,3-0; 5db: 0-5, 5-0; 10db:10-0,0-10; 15db: 15-0, 0-15 

if nargin <2
    N_c = 2;
    SNR = 2000;
end

nm1 = 'sct_esc10_16_0.25_m3_log_mix';
nm2 = ['0  ', num2str(pctrl.db), '.mat'];
nm = [nm1 nm2];    
load(nm);
data0 = data;
label0 = labels;    
nm2 = [num2str(pctrl.db), '  0.mat'];
nm = [nm1 nm2];    
load(nm);
Data = {data0, data}; % testing data with power diff
Labels = {label0, labels};


ncomb = combnk(1:10, N_c);
rng(10)
ind = randperm(200);
ind_cv = getindx(ind(1:100),featln); % get the cv data index
ind_tt = getindx(ind(101:end),featln);
cv_mixdat = [];
tt_mixdat = [];
cvmixls = [];
ttmixls = [];
for i = 1: length(Data)
    for ii = 1:size(ncomb,1) % == length(Data{i})/featln/200
    tp = Data{i}; % data
    cv_mixdat = [cv_mixdat, tp(:, (ii-1)*800 + ind_cv)];
    tt_mixdat = [tt_mixdat, tp(:, (ii-1)*800 + ind_tt)];
%     tp_l = Labels{i}; % label
%     cvmixls = [cvmixls, tp_l(:, (ii-1)*800 + ind_cv)];
%     ttmixls = [ttmixls, tp_l(:, (ii-1)*800 + ind_tt)];
    cvmixls=[cvmixls,ii*ones(1,size(ind_cv,2))];
    ttmixls=[ttmixls,ii*ones(1,size(ind_tt,2))];
    end
end









end % end of the file