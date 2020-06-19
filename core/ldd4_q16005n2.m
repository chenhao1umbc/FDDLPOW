
% loading mixture data 
if N_c == 2 % 2 class mixture
    for indx_p = 1:N_c 
        if pctrl.equal == 1 % equal power cases
            Power{indx_p} = zeros(1,N_c); % only 0 0 0 
        else
            Power{indx_p} = circshift([0, pctrl.db], indx_p-1); 
        end
    end
end

if N_c == 3 % 3 class mixture
    for indx_p = 1:N_c 
        if pctrl.equal == 1% equal power cases
            Power{indx_p} = zeros(1,N_c); % only 0 0 0 
        else
            if pctrl.if2weak == 1 %[3 0 0] the order is specific
                Power{indx_p} = circshift([pctrl.db, 0, 0], indx_p-1); 
            else
                Power{indx_p} = circshift([0, pctrl.db, pctrl.db], indx_p-1);
            end
        end
    end
end
cv_mixdat=[];tt_mixdat=[];cvmixls=[];ttmixls=[];
c = combnk(1:6, N_c); % ble bt fhss1 zb

% r00t = '.././data/RF/SNout_LMdata4/q16005n2/norm_mix';
r00t = '/extra/chenhao1/FDDLOW/data/RF/SNout_LMdata4/q16005n2/norm_mix';
% r00t = '/home/chenhao1/Matlab/FDDLOW/data/RF/SNout_LMdata4/q16005n2/norm_mix';
part1 ='db449_1classq16&0.05n2_positive_M2_snr';
part2 = 'db449_1classq16&0.05n2_negative_M2_snr';
if pctrl.if2weak == 1
part1 ='norm_db449_1classq16&0.05n2_positive_M2_snr';
part2 = 'norm_db449_1classq16&0.05n2_negative_M2_snr'; 
end
for indx_p = 1:N_c
    for indCl=1:size(c,1)
        indClnm = c(indCl, :);
        nmdb1=[r00t,num2str(indClnm), part1, num2str(SNR), 'power_',num2str(Power{indx_p}),'.mat']; 
        nmdb2=[r00t,num2str(indClnm), part2, num2str(SNR), 'power_',num2str(Power{indx_p}),'.mat'];  
        if pctrl.if2weak == 1
        nmdb1=[r00t,num2str(indClnm), part1, num2str(SNR), '.matpower_',num2str(Power{indx_p}),'.mat']; 
        nmdb2=[r00t,num2str(indClnm), part2, num2str(SNR), '.matpower_',num2str(Power{indx_p}),'.mat']; 
        end
        load(nmdb1)
        load(nmdb2)
        % concatenate the positve and negative parts
        for ii=1:length(db.features(1,:))/featln
            db2.features(:,1+(ii-1)*featln:ii*featln)=...
                flip(flip(db2.features(:,1+(ii-1)*featln:ii*featln),2),1);
        end
        db.features=[db.features;db2.features];
        % generating random indices for testing and cross-validation
        h1=1:cvln_mix;
        h2=cvln_mix+1:cvln_mix+ttln_mix;
        cvind=[]; ttind=[];
        for ii=1:featln
            cvind=[cvind;featln*h1+ii-featln];
            ttind=[ttind;featln*h2+ii-featln];
        end   
        cv_dat_temp=db.features(:,cvind); % cv samples
        tt_dat_temp=db.features(:,ttind); % testing samples
        cv_mixdat=[cv_mixdat,cv_dat_temp];
        cvmixls=[cvmixls,indCl*ones(1,size(cv_dat_temp,2))];
        tt_mixdat=[tt_mixdat,tt_dat_temp];
        ttmixls=[ttmixls,indCl*ones(1,size(tt_dat_temp,2))];
    end
end