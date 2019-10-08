function theother3 = Theother3(a,C)
% this is a funciton to get the counter-combinations
if nargin <2
    C = 6;
end
theother3 = 1:C;
for ii = 1:length(a)
    theother3(theother3 == a(ii))= [];
end

end % end of this file