close all
clear all


%discretization step
hTM = 0.1; 
%training time
Ttrain = 1000;

%initial condition for integration of Rossler system 
x0 = [-8.42496991530355,	-3.22412773526364,	0.0139822634835537];


i2 = 7; % i2    - abs() of beta degree
%reguliarization parameter
bet =  10^(-i2);

CC.lyapTime = 14.08;

plt = 0; % turn on (1)/off(0) plotting and printing intermediate results


% Total prediction time
PredictTime = 10*CC.lyapTime;   

%degree of polynomial
deg = 0:8;
%number of embedded dimension
k= 3;

% delay of embeded dimensions
tau = 1.5;


CC.tau= tau;
CC.polyn_deg = deg;
CC.k = k;

%which variable will be predicted
varbl = 2;

%Model (Rossler) parameters
a= 0.2;
b= 0.2;
c= 5.7;

P.a=a;
P.b=b;
P.c=c;


NT = ceil((Ttrain+PredictTime)/hTM)+1;
L = ceil(Ttrain/hTM)+1;


[X1,~]=Generate_signal(hTM,NT,x0,P);
T=linspace(0,hTM*(NT-1),NT);


X0 = X1(varbl,1:L);
Ntau=round(tau/hTM);
Xpred = X1(varbl,L-2*Ntau:end);
Tpred = T(L-2*Ntau:end)-T(L-2*Ntau);

[ErrorL,Terr,~,XP,TP,~,~,~,~] = Cheb_prediction(PredictTime,X0,hTM,CC,bet,Xpred,Tpred,plt);


Xmn = min(X0);
Xmx = max(X0);
a = 2/(Xmx-Xmn);
X0m = a*(Xpred-Xmn)-1;
XPm = XP;
XP = 1/a*(XP+1)+Xmn;
X0 =Xpred;

figure
subplot(311)
plot(TP/CC.lyapTime,X0), hold on
plot(TP/CC.lyapTime,XP,'--')
ylabel('$y(t)$, $\hat{y}(t)$','Interpreter','latex')
legend({'actual signal','prediction'})
xlim([0 inf])
subplot(313)
plot(TP/CC.lyapTime,abs(X0m-XPm)/2), hold on
xlabel('$\Lambda t$','Interpreter','latex')
ylabel('$\varepsilon(t)$','Interpreter','latex')
xlim([0 inf])


