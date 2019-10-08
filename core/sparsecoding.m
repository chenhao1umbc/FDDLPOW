function Z = sparsecoding(Dict, database, opts, mixture_n, cvortest)    
% this function is to control which data set to use for sparscoding
% cvortest = [1, 0] means do cv, not do test

if sum(cvortest) == 2
    error(' error from file sparsecoding.m ')
end

if mixture_n < 2

    if cvortest(1)
        Z = sparsecoding_cv(Dict, database, opts); 
    else
        Z = sparsecoding_test(Dict, database, opts);
    end
    
else  % mixture cases
    if cvortest(1)
        Z = sparsecoding_mix_cv(Dict, database, opts); 
    else
        Z = sparsecoding_mix_test(Dict, database, opts); 
    end
end


end % end of this file