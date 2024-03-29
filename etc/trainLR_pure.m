function B = trainLR_pure(opt, Database, Dict_mix)
% train logistic regression for pure cases

featln = Database.featln;
beta = opt.beta;
X = (Dict_mix.W'*Dict_mix.Z)'; % projected

lb = [1,zeros(1,5)];
for counter =1:6
    Y(featln*300*(counter-1)+1:featln*300*counter, :) = ...
    kron(ones(300*featln,1), circshift(lb, counter-1));
end

warning off
B = mnrfit(X, Y);
save(['beta', num2str(beta), 'B_X_Y_pure.mat'], 'B', 'X', 'Y');

end % end of this file
