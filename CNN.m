% {
close all
clear
clc;
tic
addpath(genpath('/extra/chenhao1/FDDLOW/fddlow'))
addpath(genpath('/extra/chenhao1/FDDLOW/data'))
addpath(genpath('/extra/chenhao1/DICTOL-master'))


% do traing or do crossvalidation
cvortest = 0;  % 1 means cv, 0 means test
SNR_INF = 2000;


% _______________ data generation __________________
% mixture_n = 2; % mixture_n classes mixture, = 1,2,3
% dbpool = [0 3 6 10 20];
% for ind = 1:5
% pctrl.db = dbpool(ind); % dynamic ratio is 0 3, 6, 10, 20 db
% if pctrl.db == 0
%     pctrl.equal = 1;
% else
%     pctrl.equal = 0;
% end
% if mixture_n < 3  pctrl.if2weak = 0; end
% if mixture_n == 3 pctrl.if2weak = 0; end
%  
% [Database]=load_data_new(mixture_n, SNR_INF, pctrl, 1000);
% Y = Database.test_mixdata;
% save(['L', num2str(mixture_n), '_', num2str(dbpool(ind)), 'db.mat'], 'Y')
% 
% end
% toc

%_____________________________get results____________________________
dbpool = [0 3 6 10 20];
mixture_n = 3; % mixture_n classes mixture, = 1,2,3
weak = zeros(1,5);
acc = zeros(1,5);
for ind = 1:5
pctrl.db = dbpool(ind); % dynamic ratio is 0 3, 6, 10, 20 db
if pctrl.db == 0
    pctrl.equal = 1;
else
    pctrl.equal = 0;
end
if mixture_n < 3  pctrl.if2weak = 0; end
if mixture_n == 3 pctrl.if2weak = 0; end
 
load(['L',num2str(mixture_n),'_', num2str(dbpool(ind)), 'db_test_ressult.mat'])
opts.C = 6; % 6 classes
opts.n = mixture_n;            

opts.equal = pctrl.equal;
if mixture_n == 3 opts.Ncombs = 20; opts.ln_test = 6000; end
if mixture_n == 2 opts.Ncombs = 15; opts.ln_test = 3000; end

[~,labels_pre] = sort(tl, 2, 'descend');
labels_pre = labels_pre';
if pctrl.if2weak == 0
    [~, weak(ind), acc(ind)] = calc_labels(labels_pre, opts);
else
    [~, weak(ind), acc(ind)] = calc_labels2w(labels_pre, opts);
end


end
toc
%%%%%%%%%%%plot results%%%%%%%%%%%%%%%%%%%

figure;
plot([0,3,6,10,20], acc, '-x')
hold on
plot([0,3,6,10,20], weak, '--o')
ylim([0.2,1])
xlabel('power ratio 0, 3,6, 10, 20db')
ylabel('accuracy')
legend('all', 'weak')
title(['L=', num2str(mixture_n), ' CNN'])



%}

