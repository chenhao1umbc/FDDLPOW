% tatal_dat is training, testing and cv pure data
dir='/home/chenhao/Matlab/LMData4/DataCollection/';
sig_name = {'ble_all', 'bt_all', 'fhss1_all', 'fhss2_all', 'wifi1_all', 'wifi2_all'};
var_name = 'x';
smp_len = 4e6; % 100 milliseconds as input lenghth
t_resol = 100;
nClass=length(sig_name);
featln= t_resol*2 -1; % time length  
fre_resol = 256;

trln=628; % trainging 300 per class
cvln=135; % cross-validation data is 50 per class
ttln=135; % testing data is 40 per class
cvttln=cvln+ttln;
trcvttln = cvttln+trln;
cvln_mix=50; % cross-validation data is 100 mixture samples per combination
ttln_mix=50; % testing data is 100 mixture samples per combination

total_dat = zeros(fre_resol, featln * (cvttln + trln) * nClass);
for ii = 1:nClass
    fname = [dir sig_name{ii}];
    obj=matfile(fname);
    len_window = floor(smp_len/t_resol);

    for k = 1:cvttln + trln
        x = obj.x(((1+(k-1)*smp_len/2):smp_len +(k-1)*smp_len/2),1);
        sp_x = spectrogram(x,len_window,[],fre_resol,4e7,'yaxis','centered');
%         temp = log(abs(sp_x));
        temp = abs(sp_x);
        total_dat(:, (1+(k-1)*featln) + featln*(cvttln+trln)*(ii-1):...
            k*featln + featln*(cvttln+trln)*(ii-1)) = temp/norm(temp, 'f');  
    end
end
save('spec_totaldat_dim256_res100','total_dat')

% mixture data part
N_c = 2;
cv_mixdat=[];tt_mixdat=[];cvmixls=[];ttmixls=[];
c = combnk(1:6, N_c); % ble bt fhss1 zb
for indx_p = 1:N_c
    for indCl=1:size(c,1)
        
        
    end
end

N_c = 3;
cv_mixdat=[];tt_mixdat=[];cvmixls=[];ttmixls=[];
c = combnk(1:6, N_c); % ble bt fhss1 zb
for indx_p = 1:N_c
    for indCl=1:size(c,1)
        
        
    end
end