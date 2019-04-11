function [acc_strong, acc_weak_av, acc_all] = calc_labels2w(labels_pre, opts)
% calculate labels for 2 weak components in 3 
% data struction is the following
% 454 356 346 ... 456 356 346 ... 456 356 346...
% section1 003    section2 300 ... section3 030
n = opts.n; % how many components in the mixture / sections
C = opts.C; % total number of classes
Ncombs = opts.Ncombs; % 6 choose 3 = 20
ln_test = opts.ln_test;

%% calculat weak component accuracy
acc_t = zeros(C, n);    
for section = 1:n 
    sect = labels_pre(:,1+(section-1)*ln_test/n:section*ln_test/n);
    % we only care about lower power accurcy
    [acc_t(:,section)] = mix_table_power(ln_test/n, section, sect, Ncombs);
end

% to calculate the weak components correct rate
function [acc0] = mix_table_power(ln_test_part, whichpart, labels_pre_part, Ncombs)            
    c = combnk(1:C,n); 
    acc0=zeros(C,1); % only score for true-positive and true negative
    for indCl=1:Ncombs
        temp = c(indCl,:);
        indClnm(1) = temp(whichpart); % find the strong component class #
        for ii0=1:ln_test_part/Ncombs          
            acc0(indClnm(1))=acc0(indClnm(1))+length(find...
                (labels_pre_part(1:n,ii0+(indCl-1)*ln_test_part/Ncombs)==indClnm(1)));            
        end
    end   
end

%% calc over-all
mixture_n = n; % sections
acc = 0;
comb = combnk(1:C, mixture_n);
ntestsample_part = ln_test/mixture_n; % n samples per section
ncomb = size(combnk(1:C, mixture_n), 1); % how many combinations 20
nsample_percomb = ntestsample_part/ncomb; % n samples percomb in a section
test_label = labels_pre(1:mixture_n, :); 
for jj = 1:mixture_n % loop over section 
    for ii = 1 + ncomb*(jj-1): ncomb*jj % loop over combinations
        for iii = 1: mixture_n % loop over each element in one combination 
            predit_labels = test_label(:, 1+nsample_percomb*(ii-1):nsample_percomb*ii);
            acc = acc + length(find(predit_labels == comb(ii-ncomb*(jj-1), iii)));
        end
    end  
end

acc_all = acc /ln_test/n;
acc_strong = sum(acc_t, 2)/nsample_percomb/size(combnk(1:(C-1),n-1),1);

if opts.equal
    acc_weak_av = acc_all;
else
    acc_weak_av = (acc_all*n-sum(acc_strong)/C)/2;
end
    
end % end of the function file