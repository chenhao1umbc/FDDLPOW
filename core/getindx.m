function [indx] = getindx(a, n_ex)
% if a == [1  3 ],  n_ex ==4, indx = [1,2,3,4, 9,10,11,12]

b = a*n_ex;
n = length(a);
indx = zeros(1,n_ex*n);
for ii = 1:n
    indx(n_ex*(ii-1)+1:ii*n_ex) = b(ii)-n_ex+1:b(ii);
end

end %end of the function file