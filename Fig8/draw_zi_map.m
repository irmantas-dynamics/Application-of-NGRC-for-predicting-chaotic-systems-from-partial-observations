close all
clear all

recPoints = [0.1; 0.5; 0.75;1;1.5]; % Errors recording points at time scale
hTM= 0.001;
Ttrain = 100; %training time

i2 = 5; % i2    - abs() of beta degree
bet =  10^(-i2);

CC.lyapTime = 1.104;

plt = 0; % turn on (1)/off(0) plotting

% Total prediction time
PredictTime = 1000.*CC.lyapTime;   

deg = 0:9;
%number of embedded dimension
k= 3;

% delay of embeded dimensions
tau = 0.1;

CC.tau= tau;
CC.k = k;
CC.polyn_deg = deg;
varbl = 3;

%Model parameters
sig=10; 
r=28;
b=8/3;
P.sig=sig; P.r=r; P.b=b;



recPointsT  =recPoints*CC.lyapTime ;

hT =hTM;
NT = ceil((Ttrain+PredictTime)/hT)+1;
L = ceil(Ttrain/hT)+1;

x0=[14.2038672459788	7.04285710095347	40.7633585468936];
[X1,~]=Generate_signal(hT,NT,x0,P);
T=linspace(0,hT*(NT-1),NT);

X0 = X1(varbl,1:L);

Ntau=round(tau/hT);
Xpred = X1(varbl,L-(k-1)*Ntau:end);
Tpred = T(L-(k-1)*Ntau:end)-T(L-(k-1)*Ntau);

[ErrorL,Terr,ErrorFixTime,~,XP,TP,~,~,~,W] = Cheb_prediction(PredictTime,recPointsT,X0,hT,CC,bet,Xpred,Tpred,plt);

M=max(max(X0),-min(X0));
Xwrite = XP'*M;

figure(1)
plot(Tpred,Xpred), hold on
plot(Tpred,Xwrite)

pid = peakfinder(Xpred',2);
pid = pid(3:end);
pks = Xpred(pid);

plot(Tpred(pid),Xpred(pid),'b*')

nn=1;
figure(2), hold on
ss=scatter(pks(1:end-nn),pks(1+nn:end),'r*');
ss.MarkerEdgeAlpha = 1;
xlabel('z_i')     
pav=sprintf('z_{i+%d}',nn);
ylabel(pav)
box on


pid = peakfinder(Xwrite',2);
pid=pid(3:end);
pks = Xwrite(pid);

figure(1),
plot(Tpred(pid),Xwrite(pid),'r*')

nn=1;
figure(2)
ss=scatter(pks(1:end-nn),pks(1+nn:end),'b.');
ss.MarkerEdgeAlpha = 1;
xlabel('z_i')     
pav=sprintf('z_{i+%d}',nn);
ylabel(pav)
box on
legend({'Signal','Prediction'});


