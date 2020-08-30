figure;
set(gcf,'color','w')

load('test_L=2_svm.mat')
p = mean(r_svm, 3);
dr = [0, 3, 6, 10, 20];
LH1 = plot(dr, nan(1, 5), '-k', dr, nan(1, 5), '--k', ...
        dr,p(1,:), '-s',dr,p(2,:),'-o',dr, p(3,:),'-^', 'linewidth', 2, 'MarkerSize', 10);
    
rw = mean(r_svm_weak, 3);

L = {'Average', 'Weak signal', 'Algorithm 1', 'Algorithm 2', 'Algorithm 3'};
legend(L)
xlabel('Power ratio (dB)')
set(gca,'FontSize',14)
set(gcf, 'Position',  [100, 100, 700, 600])
ylim([0.2,1])
ylabel('Classification accuracy')