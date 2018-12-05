clear
% close all

f_resol = 150;
t_resol = 150;
for i = 1:5
    for ii = 1:3
        figure(i*100+ii)
        filename = ['s', num2str(i),' (', num2str(ii),').wav'];
        [y, Fs] = audioread(filename);
        len_x = length(y);
        len_window = floor(len_x/t_resol);
        sp_x = spectrogram(y,len_window,[],f_resol,Fs,'yaxis','centered');
        imagesc(abs(sp_x).^2)
        title(filename)
    end
end

return

filename = ['s', num2str(1),' (', num2str(1),').wav'];
[y1, Fs] = audioread(filename);
filename = ['s', num2str(5),' (', num2str(1),').wav'];
[y5, Fs] = audioread(filename);
len_x = min( length(y1), length(y5));
len_window = floor(len_x/t_resol);
sp_x1 = spectrogram(y1(1:len_x),len_window,[],f_resol,Fs,'yaxis','centered');
sp_x5 = spectrogram(y5(1:len_x),len_window,[],f_resol,Fs,'yaxis','centered');
figure(101);imagesc(abs(sp_x1)); figure; imagesc(abs(sp_x1).^2)
figure(501);imagesc(abs(sp_x5))
% figure(601);imagesc(abs(sp_x1+sp_x1)); title('sp + sp')
figure(602);imagesc(abs(sp_x1)+abs(sp_x1)); title('abs + abs')

sp_x = spectrogram(y1(1:len_x)+y5(1:len_x),len_window,[],f_resol,Fs,'yaxis','centered');
figure(603);imagesc(abs(sp_x)); title('mix')