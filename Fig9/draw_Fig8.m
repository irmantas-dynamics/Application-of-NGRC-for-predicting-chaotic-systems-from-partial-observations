close all
clear all

recPoints = [0.1; 0.5; 0.75;1;1.5]; % Errors recording points at Lyap time scale

i2 =  9; % i2    - abs() of beta degree
bet =  10^(-i2);

CC.lyapTime = 90;

plt = 0; % turn on (1)/off(0) plotting

% Total prediction time
PredictTime = 2*CC.lyapTime;   

deg = 0:9;
%number of embedded dimension
k= 4;

hT = 1;

% delay of embeded dimensions
tau = 7*hT;


CC.tau= tau;
CC.polyn_deg = deg;
CC.k = k;

%first oscillator
varbl = 1;

recPointsT  =recPoints*CC.lyapTime ;
      

X=readmatrix('ST_0_1.dat');

X0 = X(:,varbl)';

lag = 0*1*CC.lyapTime;
Ntot = length(X0);
L = round(Ntot*0.5)+lag;
fprintf('Training signal length: %d\n',L-lag)
T = (1:Ntot)*hT;

Ntau=round(tau/hT);
LL= L-2*Ntau;

Xpred = X0(LL:end);
Tpred = T(LL:end)-T(LL);


XX0 = X0(lag+1:L);


[ErrorL,Terr,ErrorFixTime,~,XP,TP,~,~,~,~] = Cheb_prediction(PredictTime,recPointsT,XX0,hT,CC,bet,Xpred,Tpred,plt);


Xmn = min(XX0);
Xmx = max(XX0);
a = 2/(Xmx-Xmn);
X0 = a*(Xpred-Xmn)-1;


ilocMxX0 = islocalmax(X0);
ilocMxXP = islocalmax(XP);

XPResc = 1/a*(XP+1)+Xmn;

locMxX0 = Xpred(ilocMxX0);
locMxXP = XPResc(ilocMxXP);



figure
subplot(311)
plot(TP/CC.lyapTime,XPResc), hold on
plot(TP/CC.lyapTime,Xpred,'--')
ylabel('$v_2(t)$, $\hat{v}_2(t)$','Interpreter','latex')
xlim([0 5])
subplot(312)
plot(TP/CC.lyapTime,abs(X0-XP)/2), hold on
xlabel('$\Lambda t$','Interpreter','latex')
ylabel('$\varepsilon(t)$','Interpreter','latex')

xlim([0 5])

subplot(313)
plot(locMxX0(1:end-1),locMxX0(2:end),'.','MarkerSize',3,'Color',[0, 0.4470, 0.7410]), hold on
plot(locMxXP(1:end-1),locMxXP(2:end),'.','MarkerSize',3,'Color',[0.8500, 0.3250, 0.0980])

legend({'signal','prediction'})
xlabel('$v_{2,i}$','Interpreter','Latex')
ylabel('$v_{2,i+1}$','Interpreter','Latex')

% set(gcf, 'PaperSize', [8 9]);
% set(gcf,'Units','centimeters')
% set(gcf,'Position',[2,2,8,9])  



