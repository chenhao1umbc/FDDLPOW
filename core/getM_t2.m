function M = getM_t2(K, C, Nc, Z)
% This function is to get M, the mean vector per class
% H1 is block diagonal ones over how many columns per block
% H2 is all ones over how many cloumns


M=zeros(K,C);
for ii = 1:C
    M(:, ii) = mean(Z(:, 1 + Nc*(ii-1): Nc*ii), 2);
end

end % end of this function file