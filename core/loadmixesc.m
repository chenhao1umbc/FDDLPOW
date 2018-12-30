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
if pctrl.db == 0    
    nm2 = '0  0.mat';
    nm = [nm1 nm2];    
end


load(nm);







end % end of the file