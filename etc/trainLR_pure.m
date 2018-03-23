function B = trainLR_pure(opt, Database, Dict_mix)
% train logistic regression for pure cases

featln = Database.featln;
SNR = opt.SNR;

X_tr = (Dict_mix.W'*Dict_mix.Z)'; % projected

lb = [1,zeros(1,5)];
for counter =1:6
    Y_tr(featln*300*(counter-1)+1:featln*300*counter, :) = ...
    kron(ones(300*featln,1), circshift(lb, counter-1));
end

B = mnrfit(X_tr, Y_tr);
save(['SNR', num2str(SNR), 'B_X_Y_pure.mat'], 'B', 'X_tr', 'Y_tr');

end % end of this file
