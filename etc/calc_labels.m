function [acc_weak, acc_weak_av, acc_all] = calc_labels(labels_pre, opts)
    n = opts.n;
    C = opts.C;
    Ncombs = opts.Ncombs;
    ln_test = opts.ln_test;
    
    acc_t = zeros(C, n);
    % run mix_ntable
    for section = 1:n    
        sect = labels_pre(:,1+(section-1)*ln_test/n:section*ln_test/n);
        % we only care about lower power accurcy
        [acc_t(:,section)] = mix_table_power(ln_test/n, section, sect, Ncombs);
    end

    % this fucntion is shared data funcion to calc each portion of power
    function [acc0] = mix_table_power(ln_test_part, whichpart, labels_pre_part, Ncombs)            
        c = combnk(1:C,n); 
        acc0=zeros(C,1); % only score for true-positive and true negative
        for indCl=1:Ncombs
            temp = c(indCl,:);
            indClnm(1) = temp(whichpart); % find the week signal
            indClnm(2:C) = Theother3(indClnm(1)); % find non-week signal
            for ii0=1:ln_test_part/Ncombs          
                acc0(indClnm(1))=acc0(indClnm(1))+length(find...
                    (labels_pre_part(1:n,ii0+(indCl-1)*ln_test_part/Ncombs)==indClnm(1)));            
            end
        end   
    end
      
    mixture_n = n;
    acc = 0;
    comb = combnk(1:C, mixture_n);
    ntestsample_part = ln_test/mixture_n;
    ncomb = size(combnk(1:C, mixture_n), 1);
    nsample_percomb = ntestsample_part/ncomb;
    test_label = labels_pre(1:mixture_n, :); 
    for jj = 1:mixture_n
        for ii = 1 + ncomb*(jj-1): ncomb*jj
            for iii = 1: mixture_n
                predit_labels = test_label(:, 1+nsample_percomb*(ii-1):nsample_percomb*ii);
                acc = acc + length(find(predit_labels == comb(ii-ncomb*(jj-1), iii)));
            end
        end  
    end
    acc_all = acc /ln_test/n;
    acc_weak = sum(acc_t, 2)/nsample_percomb/size(combnk(1:(C-1),n-1),1);
    if opts.equal
        acc_weak_av = acc_all;
    else
        acc_weak_av = sum(acc_weak)/C;
    end
    
end % end of the function file