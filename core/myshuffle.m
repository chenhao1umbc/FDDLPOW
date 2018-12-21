function [database] = myshuffle(Old_database, seed)

% This function is made to load testing, cross-validation, and traing data with their labels
% The dataset is called GRID, we only take 10 people out of 34. 
% output is database, struct with 
%         database.featln - scattering coefficient time lenght per sample
%         database.tr_data - training data 
%         database.tr_label - training labels
%         database.cv_data - cross validation data for non-mixture
%         database.cv_label - cross validation labels
%         database.cv_mixdata - cross-val data for mixture
%         database.cv_mixlabel - cross-val labels for mixture
%         database.test_mixdata - test data for mixture
%         database.test_mixlabel -test data labels for mixture

database = Old_database;

whichclass = 1: 10;% S1 to S10
nClass = length(whichclass);
trln = 20; % trainging 300 per class
cvln = 10; % cross-validation data is 50 per class


%% loading data
featln = 4;
rng(seed)
ind = randperm(trln);
exchgind = ind(1:cvln); % get the exchange ind
exchgind_all = zeros(1, cvln * nClass * featln);
for ii = 1:nClass 
    exchgind_all(cvln*featln*(ii-1)+1:cvln*featln*ii) = getindx(exchgind,featln) + trln*featln*(ii-1); %
end
tempdata = database.cv_data;
database.cv_data = database.tr_data(:, exchgind_all);
database.tr_data(:, exchgind_all) = tempdata;

end % end of the function file