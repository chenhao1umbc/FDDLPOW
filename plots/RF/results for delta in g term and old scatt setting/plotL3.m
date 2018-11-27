weak = [0.9670	0.8703	0.6567	0.3923	0.2680 % alg 2 lr mix
0.6673	0.5870	0.4727	0.3477	0.2843 % alg 2 lr pure
0.9813	0.9267	0.7723	0.4000	0.2620 % alg 3 lr mix
0.7477	0.7083	0.6120	0.4573	0.3003 % alg 3 lr pure
0.7757	0.6990	0.6097	0.4780	0.2823 % alg 2 zf
0.8030	0.7370	0.6300	0.5083	0.2900];% alg 3 zf

all =[0.9874	0.9561	0.8854	0.7974	0.7560
0.7466	0.7421	0.7370	0.7284	0.7167
0.9907	0.9731	0.9227	0.8000	0.7540
0.7997	0.8029	0.7830	0.7530	0.7232
0.8443	0.8389	0.8276	0.8083	0.7538
0.8568	0.8539	0.8372	0.8261	0.7598];

close all
figure
for ii = 1:2
plot(weak(ii, :), ':sq')
hold on
end
for ii = 3:4
plot(weak(ii, :), ':^')
hold on
end
plot(weak(5, :), '--sq')
plot(weak(6, :), '--^')

% name = {'FDDLO with lr\_mix'
% 'FDDLO with lr\_pure'
% 'FDDLW with lr\_mix'
% 'FDDLW with lr\_pure'
% 'FDDLO with eq'
% 'FDDLW with eq'};
name = {'alg. 2 and log. regr. using mixture labels',...
    'alg. 2 and log. regr. without using mixture labels',...
    'alg. 3 and log. regr. mixture labels',...
    'alg. 3 and log. regr. without mixture labels',...
    'alg. 2 and ZF without mixture labels',...
    'alg. 3 and ZF without mixture labels'};

legend(name, 'Location','southwest')
xticklabels({'3db','', '6db','', '10db','', '20db','', '40db'})
xlabel 'Dynamic ratio in dB'
ylabel 'accuracy'
ylim([0,1])
title 'Weaker signal''s accuracy for L = 3'
set(gcf,'color','w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for ii = 1:2
plot(all(ii, :), ':sq')
hold on
end
for ii = 3:4
plot(all(ii, :), ':^')
hold on
end
plot(all(5, :), '--sq')
plot(all(6, :), '--^')

name = {'DDL with lr\_mix'
'DDL with lr\_pure'
'DDLW with lr\_mix'
'DDLW with lr\_pure'
'DDL with eq'
'DDLW with eq'};
xticklabels({'3db','', '6db','', '10db','', '20db','', '40db'})
legend(name, 'Location','southwest')
xticklabels({'3db','', '6db','', '10db','', '20db','', '40db'})
xlabel 'Dynamic ratio in dB'
ylabel 'accuracy'
ylim([0,1])
title 'accuracy for weak and non-weak signal for L = 3'

set(gcf,'color','w');