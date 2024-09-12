close all
clear all

degM = 21;

betPow = 3;

fpav = sprintf('Cheb//Cheb%d_bet%d',degM(1),betPow(1));
load(fpav)    

ErrorMCheb= ErrorM;

betPow= 9;
fpav = sprintf('Monom//Monom%d_bet%d',degM(1),betPow(1));
load(fpav)    

ErrorMMonom= ErrorM;


tt = PredictTimeM/CC.lyapTime;

medMn=zeros(length(tt),1);
medCh=zeros(length(tt),1);
q1Mn=zeros(length(tt),1);
q1Ch=zeros(length(tt),1);
q3Mn=zeros(length(tt),1);
q3Ch=zeros(length(tt),1);

for i1 = 1:length(tt)
    medMn(i1) = median(ErrorMMonom(:,i1));
    q1Mn(i1)=quantile(ErrorMMonom(:,i1),0.25);
    q3Mn(i1)=quantile(ErrorMMonom(:,i1),0.75);

    medCh(i1) = median(ErrorMCheb(:,i1));
    q1Ch(i1)=quantile(ErrorMCheb(:,i1),0.25);
    q3Ch(i1)=quantile(ErrorMCheb(:,i1),0.75);

end

figure
errorbar(tt,medCh,q1Ch-medCh,q3Ch-medCh,'LineStyle','none','Marker','_','LineWidth',1.), hold on
errorbar(tt,medMn,q1Mn-medMn,q3Mn-medMn,'LineStyle','none','Marker','_','LineWidth',1.), hold on
set(gca,'yscale','log')

xlim([0 0.95])
ylim([2e-4 0.3])
xticks(0.1:0.2:0.9)
ylabel('RSME','Interpreter','latex')
xlabel('$T_{\mathrm{pred}}\Lambda$','Interpreter','latex')
legend({'Ch21','Mn21'})

set(gcf, 'PaperSize', [9 7]);
set(gcf,'Units','centimeters')
set(gcf,'Position',[2,2,9,7])  



