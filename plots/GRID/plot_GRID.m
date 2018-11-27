clear

f_resol = 500;
t_resol = 10;
for i = 1:3
    for ii = 1:5
        figure(i*100+ii)
        filename = ['s', num2str(i),' (', num2str(ii),').wav'];
        [y, Fs] = audioread(filename);
        len_x = length(y);
        len_window = floor(len_x/t_resol);
        sp_x = spectrogram(y,len_window,[],f_resol,Fs,'yaxis','centered');
        imagesc(abs(sp_x))
        title(filename)
    end
end



