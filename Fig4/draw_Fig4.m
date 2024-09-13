close all
clear all
load("results2.mat")

nn=1;
figure
subplot(211)
plot(XsignalMaxArray(1:nn:end,1),XsignalMaxArray(1:nn:end,2),'.','Color',[0 0.4470 0.7410],'MarkerSize',1), 
hold on
for i1 =1:length(pM)
    plot([1 1]*pM(i1),[3 7],'--','Color',[0.4660 0.6740 0.1880])
end
xlim([0.58 0.75])
ylim([3 7])

ylabel('$y$','Interpreter','latex');
box on

subplot(212)
plot(XpredMaxArray(1:nn:end,1),XpredMaxArray(1:nn:end,2),'.','Color',[0.8500 0.3250 0.0980],'MarkerSize',1)
hold on
box on
for i1 =1:length(pM)
    plot([1 1]*pM(i1),[3 7],'--','Color',[0.4660 0.6740 0.1880])
end
ylim([3 7])

xlim([0.58 0.75])
xlabel('$b$','Interpreter','latex');
ylabel('$\hat{y}$','Interpreter','latex');
set(gcf, 'PaperSize', [8 7]);
set(gcf,'Units','centimeters')
set(gcf,'Position',[2,2,8,7])  