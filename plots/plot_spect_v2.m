% close all
clear

smp_len = 4e6;
k = 3; % which sample

dir='/home/chenhao/Matlab/LMData4/DataCollection/';
sig_name = {'ble_all', 'bt_all', 'fhss1_all', 'fhss2_all', 'wifi1_all', 'wifi2_all'};
var_name = 'x';
whichsignal = 2;
t_resol = 100;

fname = [dir sig_name{whichsignal}];
obj=matfile(fname);
x=obj.(var_name)(((1+(k-1)*smp_len):k*smp_len),1);
% x=circshift(x,1e6);

len_x = length(x); % xlen
len_window = floor(len_x/t_resol);
%nov = floor(nsc/2); % hop
%nff = max(256,2^nextpow2(nsc));% nfft number of fft points 

% figure
% spectrogram(x,len_window,[],1024,4e7,'yaxis','centered'); % 4e7 is the sampling frequency
sp_x = spectrogram(x,len_window,[],256,4e7,'yaxis','centered');

figure
% imagesc(log(abs(sp_x)))
imagesc(abs(sp_x))
title(sig_name{whichsignal})