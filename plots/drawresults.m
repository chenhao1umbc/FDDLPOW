close all

%{
FDDL=[0.9993
0.99166
0.9965
0.969833333
0.816833333];

FDDLO=[
1
1
1
0.9957
0.8974];

FDDLOW = [
1
1
0.9998
0.9948
0.8828
];

L2_old_all = [0.9955 0.9780 	0.9860 	0.9867 	0.9347 ];

figure
plot(FDDL,'x-')
hold on
plot(FDDLO,'x-')
plot(FDDLOW, 'x-')
plot(L2_old_all, 'x-')

legend('FDDL','FDDLO', 'FDDLOW','old\_FDDLOW')
grid on
grid minor
ylim([0.4, 1.02])
title ('Overall signal accuracy for L = 2')
xticklabels({'0db','','3db','','6db','','10db','','20db'})
%}

%%%%%%%%%%% L = 3%%%%%%%%%%%%%%
FDDL3=[0.9657
0.9592
0.9396
0.8984
0.8142];

FDDLO3=[
0.9813
0.9797
0.9711
0.9445
0.8386
];

FDDLOW3 = [
0.9812
0.9797
0.9711
0.9443
0.8393];

% L3_all_old = [ 0.8755 0.8568	0.8539	0.8372	0.8261];

figure
plot(FDDL3,'x-')
hold on
plot(FDDLO3,'x-')
plot(FDDLOW3, 'x-')
% plot(L3_all_old, 'x-')

legend('FDDL','FDDLO', 'FDDLOW')
grid on
grid minor
ylim([0.4, 1.02])
title ('Overall signal accuracy for L = 3')
xticklabels({'0db','','3db','','6db','','10db','','20db'})

%% week signal
%{
FDDL=[
0.9993
0.9983
0.993
0.9397
0.6337];


FDDLO=[
1
1
1
0.9913
0.7893];

FDDLOW = [
1
1
0.9997
0.9897
0.7657];


L2_old = [0.97 0.9613 	0.9567 	0.9493 	0.8987 ];

figure
plot(FDDL,'x-')
hold on
plot(FDDLO,'x-')
plot(FDDLOW, 'x-')
plot(L2_old, 'x-')

legend('FDDL','FDDLO', 'FDDLOW', 'old\_FDDLOW')
grid on
grid minor
ylim([0.4, 1.02])
title ('Week signal accuracy for L = 2')
xticklabels({'0db','','3db','','6db','','10db','','20db'})
%}


FDDL3=[0.9657
0.927
0.848
0.6988
0.4468
];

FDDLO3=[
0.9813
0.9732
0.9348
0.8335
0.5158];

FDDLOW3 = [
0.9812
0.9723
0.9355
0.833
0.5178];


% L3_old = [0.85  0.8030	0.7370	0.6300	0.5083];

figure
plot(FDDL3,'x-')
hold on
plot(FDDLO3,'x-')
plot(FDDLOW3, 'x-')
% plot(L3_old, 'x-')

legend('FDDL','FDDLO', 'FDDLOW')
grid on
grid minor
ylim([0.4, 1.02])
title ('Week signal accuracy for L = 3')
xticklabels({'0db','','3db','','6db','','10db','','20db'})

%% old vs new
%{
L3_old = [0.85  0.8030	0.7370	0.6300	0.5083];
L3 = [0.9813    0.9732    0.9348    0.8335    0.5158];

figure
plot(L3,'x-')
hold on
plot(L3_old,'x-')

ylim([0, 1.05])
legend('new','old')
title ('Week signal for L = 3')
xticklabels({'0db','','3db','','6db','','10db','','20db'})
grid minor
grid on


L2_old = [0.97 0.9613 	0.9567 	0.9493 	0.8987 ];
L2 = [1 1 1 0.9913 0.7893 ];

figure
plot(L2,'x-')
hold on
plot(L2_old,'x-')

ylim([0, 1.05])
legend('new','old')
title ('Week signal for L = 2')
xticklabels({'0db','','3db','','6db','','10db','','20db'})
grid minor
grid on
%}