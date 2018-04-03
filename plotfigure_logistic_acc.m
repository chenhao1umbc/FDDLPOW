

lgwithmix_alg2 = [0.9700 	0.9527 	0.8127 	0.3353 	0.2000 ];
lgwithmix_alg3 = [0.9980 	0.9960 	0.9560 	0.4180 	0.2047 ];
lgnomix_alg2 = [0.8440 	0.7993 	0.7020 	0.4660 	0.2187 ];
lgnomix_alg3 = [0.9827 	0.9773 	0.8947 	0.5380 	0.2767 ];
eqnomix_alg2 = [0.9580 	0.9580 	0.9600 	0.9007 	0.5707 ];
eqnomix_alg3 = [0.9400 	0.9360 	0.9187 	0.8587 	0.5473 ];
db = 1:5;
figure 
plot(db,lgwithmix_alg2, '--sq')
xticklabels({'3', '', '6', '', '10','', '20','', '40'})
hold on
plot(db,lgwithmix_alg3, '--^')
plot(db,lgnomix_alg2, '--sq')
plot(db,lgnomix_alg3, '--^')
plot(db,eqnomix_alg2, ':sq')
plot(db,eqnomix_alg3, ':^')
legend('alg. 2 and log. regr. using mixture labels',...
    'alg. 3 and log. regr. using mixture labels',...
    'alg. 2 and log. regr. without mixture labels',...
    'alg. 3 and log. regr. without mixture labels',...
    'alg. 2 and ZF without mixture labels',...
    'alg. 3 and ZF without mixture labels','Location', 'southwest');
title 'Weaker signal accuracy'
xlabel('Dynamic ratio in dB')
ylabel 'accuracy for weaker signal'
set(gcf,'color','w');
