function [normdata] = norm_data(data)
% normalize the data per column

normdata = data./sqrt(sum(data.^2,1));
normdata(isnan(normdata)) = 0;
end %end of this file