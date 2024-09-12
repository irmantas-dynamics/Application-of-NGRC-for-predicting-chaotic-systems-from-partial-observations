% close all
clear all


recPoints = [0.1; 0.5; 0.75;1;1.5]; % Errors recording points at time scale
hTM= 0.005;
Ttrain = 100; %training time


i2 = 9; % i2    - abs() of beta degree
bet =  10^(-i2);

CC.lyapTime = 1.104;

plt = 0; % turn on (1)/off(0) plotting


% Total prediction time
PredictTime = 3*CC.lyapTime;   

deg = 1:2:9;
%number of embedded dimension
k= 3;

% delay of embeded dimensions
tau = 0.15;

CC.tau= tau;
CC.k = k;
CC.polyn_deg = deg;
varbl = 1;

%Model parameters
sig=10; 
r=28;
b=8/3;
P.sig=sig; P.r=r; P.b=b;



recPointsT  =recPoints*CC.lyapTime ;

hT =hTM;
NT = ceil((Ttrain+PredictTime)/hT)+1;
L = ceil(Ttrain/hT)+1;


x0= [14.2038672459788	7.04285710095347	40.7633585468936];
[X1,~]=Generate_signal(hT,NT,x0,P);
 T=linspace(0,hT*(NT-1),NT);


X0 = X1(varbl,1:L);
Ntau=round(tau/hT);
Xpred = X1(varbl,L-2*Ntau:end);
Tpred = T(L-2*Ntau:end)-T(L-2*Ntau);



[ErrorL,Terr,ErrorFixTime,~,XP,TP,~,~,~,~] = Monom_prediction(PredictTime,recPointsT,X0,hT,CC,bet,Xpred,Tpred,plt);


M=max(max(X0),-min(X0));
X0 =Xpred/M;

subplot(312)
plot(TP/CC.lyapTime,X0*M), hold on
plot(TP/CC.lyapTime,XP*M,'--')
ylabel('$x(t)$, $\hat{x}(t)$','Interpreter','latex')
xlim([0 inf])
subplot(313)
plot(TP/CC.lyapTime,abs(X0-XP)/2), hold on
xlabel('$\Lambda_1 t$','Interpreter','latex')
ylabel('$\varepsilon(t)$','Interpreter','latex')
xlim([0 3])
legend({'Ch9','Mn9'})


