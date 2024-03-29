function Z = sbrvm(Dict, Database, mixture_n, cvortest)
% this function is to control which data set to use for sparscoding
% cvortest = [1, 0] means do cv, not do test


if cvortest == 1 
    if mixture_n == 1
        data = Database.cv_data;
    else
        data = Database.cv_mixdata;
    end
    
else % validation cases
    if mixture_n == 1
        data = Database.test_data;
    else
        data = Database.test_mixdata;
    end
end

n = size(data, 2);
Z = zeros(size(Dict.D, 2), n);
p = cell(1, n);
parfor i = 1:n
    [p{i}, ~, ~] = SparseBayes('Gaussian', Dict.D, data(:,i));
end
for i = 1:n
Z(p{i}.Relevant, i) = p{i}.Value;
end



end % end of this file