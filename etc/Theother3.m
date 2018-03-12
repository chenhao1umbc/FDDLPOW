function theother3 = Theother3(a)
% this is a funciton to get the counter-combinations
theother3 = 1:6;
for ii = 1:length(a)
    theother3(theother3 == a(ii))= [];
end

end % end of this file