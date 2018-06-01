function [M, H1, H2] = getMH1H2_t2(trlabels, Z)
% This function is to get M, the mean vector per class
% H1 is block diagonal ones over how many columns per block
% H2 is all ones over how many cloumns

C = max(trlabels);
[K, N]=size(Z);
Nc = N / C;
M=zeros(K,C);
for ii = 1:C
    M(:, ii) = mean(Z(:, 1 + Nc*(ii-1): Nc*ii), 2);
end
H1 = kron(eye(C),ones(Nc)/Nc);
H2 = ones(N)/N;

end % end of this function file