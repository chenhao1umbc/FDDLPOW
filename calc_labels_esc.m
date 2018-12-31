function [acc] = calc_labels_esc(labels_pre, opts)
    mixture_n = opts.n;    
    acc = 0;
    comb = combnk(1:opts.C, mixture_n);
    ntestsample_part = opts.ln_test/mixture_n;
    ncomb = size(combnk(1:opts.C, mixture_n), 1);
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
    acc_all = acc /opts.ln_test/n;
    
end % end of the function file