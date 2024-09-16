close all
clear all

%maximum degree of polynomial (in publication marked as m)
degM = 14;

%degree of regression parameter, i.e. beta = 10^(-betaPow)
betPow = 9;

%loading data of Chebyshev polynomials
fpav = sprintf('Cheb//Cheb%d_bet%d',degM(1),betPow(1));
load(fpav)    
ErrorMCheb= ErrorM;

%loading data of polynomials
betPow= 9;
fpav = sprintf('Monom//Monom%d_bet%d',degM(1),betPow(1));
load(fpav)    
ErrorMMonom= ErrorM;

%array of continuous prediction time
tt = PredictTimeM/CC.lyapTime;

medMn=zeros(length(tt),1);
medCh=zeros(length(tt),1);
q1Mn=zeros(length(tt),1);
q1Ch=zeros(length(tt),1);
q3Mn=zeros(length(tt),1);
q3Ch=zeros(length(tt),1);


for i1 = 1:length(tt)
    %median, first and third quartiles of continuous monomial prediction
    %error
    medMn(i1) = median(ErrorMMonom(:,i1));
    q1Mn(i1)=quantile(ErrorMMonom(:,i1),0.25);
    q3Mn(i1)=quantile(ErrorMMonom(:,i1),0.75);

    %median, first and third quartiles of continuous Chebyshev polynomial prediction
    %error
    medCh(i1) = median(ErrorMCheb(:,i1));
    q1Ch(i1)=quantile(ErrorMCheb(:,i1),0.25);
    q3Ch(i1)=quantile(ErrorMCheb(:,i1),0.75);

end

figure
errorbar(tt,medCh,q1Ch-medCh,q3Ch-medCh,'LineStyle','none','Marker','_','LineWidth',1.), hold on
errorbar(tt,medMn,q1Mn-medMn,q3Mn-medMn,'LineStyle','none','Marker','_','LineWidth',1.), hold on

xlim([1 9])
ylabel('RSME','Interpreter','latex')
set(gca,'yscale','log')
xlabel('$T_{\mathrm{pred}}\Lambda$','Interpreter','latex')
legend({'Ch14','Mn14'})

set(gcf, 'PaperSize', [9 7]);
set(gcf,'Units','centimeters')
set(gcf,'Position',[2,2,9,7])  



