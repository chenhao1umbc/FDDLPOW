close all
clear

smp_len = 8e6;
k = 3; % which sample

dir='/home/chenhao/Matlab/LMData4/DataCollection/';
sig_name = {'ble_all', 'bt_all', 'fhss1_all', 'fhss2_all', 'wifi1_all', 'wifi2_all'};
var_name = 'x';
whichsignal = 1;
f_resol = 256;
t_resol = 50;

% signal 1
fname = [dir sig_name{whichsignal}];
obj=matfile(fname);
x=obj.(var_name)(((1+(k-1)*smp_len):k*smp_len),1);
x = x/norm(x);

% signal 2
fname = [dir sig_name{2}];
obj=matfile(fname);
x2=obj.(var_name)(((1+(k-1)*smp_len):k*smp_len),1);
x2 = x2/norm(x2);

len_x = length(x); % xlen
len_window = floor(len_x/t_resol);
%nov = floor(nsc/2); % hop
%nff = max(256,2^nextpow2(nsc));% nfft number of fft points 


%% plot figures

sp_x = spectrogram(x,len_window,[],f_resol,4e7,'yaxis','centered');
sp_x2 = spectrogram(x2,len_window,[],f_resol,4e7,'yaxis','centered');
sp_xplus = spectrogram(x+x2,len_window,[],f_resol,4e7,'yaxis','centered');

figure
% imagesc(log(abs(sp_x)))
% title(['singal\_1', ' log'])
imagesc(abs(sp_x))
title(['singal\_1', ' without\_log'])
% caxis([-18, 1])
colorbar
title(['200ms\_BLE : ', num2str(f_resol), ' freq. & ', num2str(2*t_resol-1),' time'])
%{
figure
% imagesc(log(abs(sp_x2)))
% title(['signal\_2', ' log'])
imagesc(abs(sp_x2))
title(['signal\_2', ' without\_log'])
% caxis([-18, 1])
colorbar

figure
% imagesc(log(abs(sp_xplus)))
imagesc(abs(sp_xplus))
title([' sum before spectrogram'])
% caxis([-18, 1])
colorbar

figure
% imagesc(log(abs(sp_x2))+log(abs(sp_x)))
imagesc(abs(sp_x2)+abs(sp_x))
title([' sum of spectrogram'])
% caxis([-18, 1])
colorbar
%}