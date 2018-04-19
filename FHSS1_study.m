
% all the codes show projection for L =2
symb = {'x','sq', '*', '>','^','.'};
for ii = 1:6
rt(ii, :) = H(:,ii)'*W'*aoos(Z,4,size(Z, 2))/norm(H(:, ii),2);
% figure; 
plot(rt(ii, 1:end/2), symb{ii}); grid on; hold on
end
legend('ble', 'bt', 'fhss1', 'fhss2', 'wif1', 'wifi2')

pr = rt(:, [251:300,1001:1050]);
[~, lb] = sort(pr, 1, 'descend');
imagesc(pr)
figure
imagesc(Dict_mix.D); title 'Dictionary'; colorbar


Zs = aoos(Z,4,size(Z, 2));
tr_Zs = aoos(Dict_mix.Z,4,size(Dict_mix.Z, 2));
figure; imagesc([W'*Zs(:, [251:300,1001:1050]), W'*tr_Zs(:, [301:900, 1501:1800])])

Z2 = tr_Zs(:, 301:600);
Z3 = tr_Zs(:, 601:900);
Z6 = tr_Zs(:, 1501:1800);
figure; imagesc(Z2+ Z3); title 'trained BT+FHSS1';colorbar
figure; imagesc(Z2+ Z6); title 'trained BT+Wifi2';colorbar

Z_weak2_3 = Zs(:, 251:300); % equal power case Z_weak2_3 == Z_2_weak3
Z_2_weak3 = Zs(:, 1001:1050);

figure; imagesc(Z2); title 'bt'; colorbar
figure; imagesc(Z3); title 'FHSS1'; colorbar
figure; imagesc(Z6); title 'wifi2'; colorbar
figure; imagesc(Z_weak2_3); title 'test bt + FHSS1'; colorbar
% figure; imagesc(Z_weak2_3); title 'weak bt + FHSS1'; colorbar
% figure; imagesc(Z_2_weak3); title 'bt + weak FHSS1'; colorbar
close all

figure; imagesc(W'*Z2); title 'bt'; colorbar
figure; imagesc(W'*Z3); title 'FHSS1'; colorbar
figure; imagesc(W'*Z6); title 'wifi2'; colorbar
figure; imagesc(W'*Z_weak2_3); title 'test bt + FHSS1'; colorbar
% figure; imagesc(W'*Z_weak2_3); title 'weak bt + FHSS1'; colorbar
% figure; imagesc(W'*Z_2_weak3); title 'bt + weak FHSS1'; colorbar
figure; imagesc(W'*(Z2 + Z3)); title 'bt + FHSS1'; colorbar
figure; imagesc(W'*(Z2 + Z6)); title 'bt + wifi2'; colorbar
close all

Z0 = (3.16*Z2+Z3);
result0 = pinv(H)*W'*aoos(Z0,4,size(Z0, 2));
[~, labels_pre0] = sort(result0, 1, 'descend');

Z0 = (10* Z2+Z3);
result0 = pinv(H)*W'*aoos(Z0,4,size(Z0, 2));
[~, labels_pre0] = sort(result0, 1, 'descend');
