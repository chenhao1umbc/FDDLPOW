

lgwithmix_alg2 = [0.9700 	0.9527 	0.8127 	0.3353 	0.2000 ];
lgwithmix_alg3 = [0.9980 	0.9960 	0.9560 	0.4180 	0.2047 ];
lgnomix_alg2 = [0.8440 	0.7993 	0.7020 	0.4660 	0.2187 ];
lgnomix_alg3 = [0.9827 	0.9773 	0.8947 	0.5380 	0.2767 ];
db = 1:5;
figure 
plot(db,lgwithmix_alg2, '--sq')
xticks([3, 6, 10, 20, 40])
hold on
plot(db,lgwithmix_alg3, '--^')
plot(db,lgnomix_alg2, '--sq')
plot(db,lgnomix_alg3, '--^')
legend('algorithem 2 using mixture labels',...
    'algorithem 3 using mixture labels',...
    'algorithem 2 without using mixture labels',...
    'algortithem 3 without using mixture labels')
title 'algorithem 2 vs 3 with&without mixture label information'