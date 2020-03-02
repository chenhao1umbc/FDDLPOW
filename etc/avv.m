
function ar_zf = avv(r_zf, ar_zf)
    for mixture_n = 1:3
    for indd = [1,2]   
    for alg = 1:3
        s = 0;
        for f = 1:5
            if ~isempty(r_zf{mixture_n, indd, f, alg})
                s = s + r_zf{mixture_n, indd, f, alg};    
            end
        end
        ar_zf{mixture_n, indd, alg} = s/5;
    end
    end
    end
end % end of the function