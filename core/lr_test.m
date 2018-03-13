function [acc, acc_av] = lr_test(Dict, Database, Z, B)

C = 6; % 6 classes
featln = 4;
n = Database.N_c;
N = size(Database.test_mixlabel, 2);
ln_test= N/featln;
W = Dict.W;
wz = (W'*aoos(Z, featln, N));
pre_prob = mnrval(B, wz');
[~,labels_pre] = sort(pre_prob, 2, 'descend');
labels_pre = labels_pre';

acc_t = zeros(C, n);
% run mix_ntable
for section = 1:n    
    sect = labels_pre(:,1+(section-1)*ln_test/n:section*ln_test/n);
    % we only care about lower power accurcy
    [acc_t(:,section)] = mix_table_power(ln_test/n, section, sect);
end

acc = sum(acc_t, 2)/50/size(combnk(1:5,n-1),1);
acc_av = sum(acc)/C;

% this fucntion is shared data funcion to calc each portion of power
function [acc0] = mix_table_power(ln_test_part, whichpart,labels_pre_part)        
    Ncombs=max(Database.cv_mixlabel);
    c = combnk(1:6,n); % ble bt fhss1 zb
    acc0=zeros(C,1); % only score for true-positive and true negative

    for indCl=1:Ncombs
        temp = c(indCl,:);
        indClnm(1) = temp(whichpart);
        indClnm(2:6) = Theother3(indClnm(1));
        for ii=1:ln_test_part/Ncombs          
                acc0(indClnm(1))=acc0(indClnm(1))+length(find...
                    (labels_pre_part(1:n,ii+(indCl-1)*ln_test_part/Ncombs)==indClnm(1)));            
        end

    end   
   
end


end % end of the function file
