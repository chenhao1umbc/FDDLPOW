function plot_rr_(r1, r1_weak)
rr = r1;
for a = 1:3
for d = 1:5
p(a, d) = (sum(rr(a, d, :)) - max(rr(a, d, :)) -min(rr(a, d, :)))/3;
end
end
p
% p = mean(r1, 3);

dr = [0, 3, 6, 10, 20];
figure;
set(gcf,'color','w')
LH1 = plot(dr, nan(1, 5), '-k', dr, nan(1, 5), '--k', ...
        dr,p(1,:), '-s',dr,p(2,:),'-o',dr, p(3,:),'-^', 'linewidth', 2, 'MarkerSize', 10);

rr = r1_weak;
for a = 1:3
for d = 1:5
p(a, d) = (sum(rr(a, d, :)) - max(rr(a, d, :)) -min(rr(a, d, :)))/3;
end
end
p
% p = mean(rr, 3);

hold on
LH3 = plot(dr,p(1,:), '--sb',dr,p(2,:),'--or',dr, p(3,:),'--^g','linewidth', 2, 'MarkerSize', 10);
ylim([0.2,1])
xlabel('Power ratio (dB)')
set(gca,'FontSize',14)
set(gcf, 'Position',  [100, 100, 700, 600])
ylabel('Classification accuracy')


LH(1) = plot(nan, nan, 'k-');
LH(2) = plot(nan, nan, 'k--');
L = {'Average', 'Weak signal', 'Algorithm 1', 'Algorithm 2', 'Algorithm 3'};


legend(LH1, L)
% legend({'Algorithm 1', 'Algorithm 2', 'Algorithm 3', 'Algorithm 1 - weak', 'Algorithm 2 - weak', 'Algorithm 3 - weak'})
end