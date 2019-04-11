
close all
clear
clc;

addpath(genpath('.././fddlow'))
addpath(genpath('.././data'))

%% Results Plot L = 2 for alg1 alg2 alg3
FDDL=[0.9993 0.99166 0.9965 0.969833333 0.816833333];
FDDLO=[1 1 1 0.9957 0.8974];
FDDLOW = [1 1 0.9998 0.9948 0.9128];
% week signal
FDDL_=[0.9993 0.9983 0.993 0.9397 0.6337];
FDDLO_=[1 1 1 0.9913 0.7893];
FDDLOW_= [1 1 0.9997 0.9897 0.8152];

figure
plot(FDDL,'rx-')
hold on
plot(FDDLO,'gx-')
plot(FDDLOW, 'bx-')

plot(FDDL_,'r^--')
plot(FDDLO_,'g^--')
plot(FDDLOW_, 'b^--')

legend('alg. 1 and Zero-forcing classifier, overall accuracy',...
    'alg. 2 and Zero-forcing classifier, overall accuracy',...
    'alg. 3 and Zero-forcing classifier, overall accuracy',...
    'alg. 1 and Zero-forcing classifier for weak component',...
    'alg. 2 and Zero-forcing classifier for weak component',...
    'alg. 3 and Zero-forcing classifier for weak component','Location', 'southwest');
grid on
% grid minor
ylim([0.4, 1.02])
title ('L = 2')
xticklabels({'0db','','3db','','6db','','10db','','20db'})
xlabel('Dynamic ratio in dB')
ylabel 'Accuracy'
set(gcf,'color','w');

%% Results Plot L = 3 with 2 weak components for alg1 alg2 alg3 
alg1 = [0.9651, 0.9600, 0.9443, 0.9092, 0.8156];
alg1_w = [0.9651, 0.9493, 0.9164,0.8638, 0.7233];

alg2 = [0.9816, 0.9812, 0.9802, 0.9683, 0.8992];
alg2_w = [0.9816, 0.9784, 0.9727,0.9425, 0.8288];

alg3 = [0.9812, 0.9793, 0.9779, 0.9684, 0.9050];
alg3_w = [0.9812, 0.9773, 0.9713,0.9526, 0.8575];


figure
plot(alg1,'rx-')
hold on
plot(alg2,'gx-')
plot(alg3, 'bx-')

plot(alg1_w,'r^--')
plot(alg2_w,'g^--')
plot(alg3_w, 'b^--')

legend('alg. 1 and Zero-forcing classifier, overall accuracy',...
    'alg. 2 and Zero-forcing classifier, overall accuracy',...
    'alg. 3 and Zero-forcing classifier, overall accuracy',...
    'alg. 1 and Zero-forcing classifier for weak component',...
    'alg. 2 and Zero-forcing classifier for weak component',...
    'alg. 3 and Zero-forcing classifier for weak component','Location', 'southwest');
grid on
% grid minor
% ylim([0.4, 1.02])
title ('L = 3 with two weak components')
xticklabels({'0db','','3db','','6db','','10db','','20db'})
xlabel('Dynamic ratio in dB')
ylabel 'Accuracy'
set(gcf,'color','w');

%% %% Results Plot L = 3 for alg1 alg2 alg3
alg1 = [0.9651, 0.9579, 0.9386, 0.8988, 0.8143];
alg1_w = [0.9651, 0.9233, 0.8440,0.6988, 0.4428];

alg2 = [0.9816, 0.9797, 0.9710, 0.9369, 0.8401];
alg2_w = [0.9816, 0.9723, 0.9350,0.8155, 0.5103];

alg3 = [0.9812, 0.9794, 0.9633, 0.9434, 0.8446];
alg3_w = [0.9812, 0.9720, 0.9215,0.8303, 0.5337];


figure
plot(alg1,'rx-')
hold on
plot(alg2,'gx-')
plot(alg3, 'bx-')

plot(alg1_w,'r^--')
plot(alg2_w,'g^--')
plot(alg3_w, 'b^--')

legend('alg. 1 and Zero-forcing classifier, overall accuracy',...
    'alg. 2 and Zero-forcing classifier, overall accuracy',...
    'alg. 3 and Zero-forcing classifier, overall accuracy',...
    'alg. 1 and Zero-forcing classifier for weak component',...
    'alg. 2 and Zero-forcing classifier for weak component',...
    'alg. 3 and Zero-forcing classifier for weak component','Location', 'southwest');
grid on
% grid minor
% ylim([0.4, 1.02])
title ('L = 3')
xticklabels({'0db','','3db','','6db','','10db','','20db'})
xlabel('Dynamic ratio in dB')
ylabel 'Accuracy'
set(gcf,'color','w');




%% Plot the training samples for each class
%{
[Database]=load_data_new(1, 2000);
name = {'BLE','BT','FHSS1','FHSS2','WIFI1','WIFI2'};
for i = 1:6
figure
imagesc(Database.tr_data(:,1+14400/6*(i-1):14400/6*(i-1)+800))
title(['Deep scattering spectrum of ', name{i},' training data'])
xlabel('Sample index')
ylabel('Frequency index')
set(gcf,'color','w');
end
%}

%% Using the logistic regression vs Zero-Forcing
%{
lgwithmix_alg2 = [0.9700 	0.9527 	0.8127 	0.3353 	0.2000 ];
lgwithmix_alg3 = [0.9980 	0.9960 	0.9560 	0.4180 	0.2047 ];
lgnomix_alg2 = [0.8440 	0.7993 	0.7020 	0.4660 	0.2187 ];
lgnomix_alg3 = [0.9827 	0.9773 	0.8947 	0.5380 	0.2767 ];
eqnomix_alg2 = [0.9580 	0.9580 	0.9600 	0.9007 	0.5707 ];
% eqnomix_alg3 = [0.9400 	0.9360 	0.9187 	0.8587 	0.5473 ];
eqnomix_alg3 =[0.9600 	0.9587 	0.9607 	0.8927 	0.5280];

db = 1:5;
figure 
plot(db,lgwithmix_alg2, ':sq')
xticklabels({'3', '', '6', '', '10','', '20','', '40'})
hold on
plot(db,lgnomix_alg2, ':sq')
plot(db,lgwithmix_alg3, ':^')
plot(db,lgnomix_alg3, ':^')
plot(db,eqnomix_alg2, '--sq')
plot(db,eqnomix_alg3, '--^')
legend('alg. 2 and log. regr. using mixture labels',...
    'alg. 2 and log. regr. without using mixture labels',...
    'alg. 3 and log. regr. mixture labels',...
    'alg. 3 and log. regr. without mixture labels',...
    'alg. 2 and ZF without mixture labels',...
    'alg. 3 and ZF without mixture labels','Location', 'southwest');

title 'Weaker signal''s accuracy  L = 2'
xlabel('Dynamic ratio in dB')
ylabel 'accuracy for weaker signal'
ylim([0,1])
set(gcf,'color','w');
%}