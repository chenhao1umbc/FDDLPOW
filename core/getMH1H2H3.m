function [H1, H2, H3] = getMH1H2H3(N, Nc, C)
% This function is to get M, the mean vector per class
% H1 is block diagonal ones over how many columns per block
% H2 is all ones over how many cloumns
% H3 is defined from M = Z*H3, M is the mean of each class

H1 = kron(eye(C),ones(Nc)/Nc);
H2 = ones(N)/N;
H3 = kron(eye(C),ones(Nc, 1)/Nc);

end % end of this function file