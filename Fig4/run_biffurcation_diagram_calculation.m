close all
clear all

hT = 0.1;
Ttrain = 1000; %training time

i2 = 3; % i2    - abs() of beta degree
bet =  10^(-i2);

CC.lyapTime = 1;

plt = 0; % turn on (1)/off(0) plotting



% Total prediction time
PredictTime = 100*CC.lyapTime;   

deg = 0:7;
%number of embedded dimension
k= 4;

% delay of embeded dimensions
tau = 1.5;


CC.tau= tau;
CC.polyn_deg = deg;
CC.k = k;

%varbl = [1,2];
varbl = 2;
CC.dim = length(varbl);
CC.dim_c = 1;
CC.Cmn =0;
CC.Cmx = 1;
CC.varbl = varbl;

%Model parameters
a= 0.2;
b= 0.4; % nestandartinis parametras
c= 5.7;

P.a=a;
P.b=b;
P.c=c;

h= hT;
dim = CC.dim;
dim_c =CC.dim_c;

recPointsT = [0.1; 0.5; 0.75;1;1.5]; % Errors recording
    
NT = ceil((Ttrain+PredictTime)/hT)+1;
L = ceil(Ttrain/hT)+1;


x0=[-3.07972019640287	-0.410252884196477	0.0228514537795281];

pM = linspace(0.6,0.65,4);

[Xtrain,Dtrain,resc,Xmn] = prepare_train_signal(hT,L,x0,P,pM,varbl);

% plotting training data    
%     figure
%     plot(Xtrain(1,:)), hold on
%     for i1 =1:length(pM)
%         plot([i1*L i1*L],[-1 1],'r-')
%     end

P.b= b;

[X1,~]=Generate_signal(hT,NT,x0,P);
x00=X1(:,end);

[X1,~]=Generate_signal(hT,NT,x00,P);
x00=x0;
T=linspace(0,hT*(NT-1),NT);

Ntau=round(tau/hT);
LL = length(X1(1,L-2*Ntau:end));
Xpred = zeros(dim+dim_c,LL);
for i1 = 1:dim
    tmpX = X1(varbl(i1),L-2*Ntau:end);
    Xpred(i1,:) = resc(i1)*(tmpX-Xmn(i1))-1;
end
Xpred(dim+1,:) = ones(1,LL)*P.b;
Tpred = T(L-2*Ntau:end)-T(L-2*Ntau);

[ErrorL,Terr,ErrorFixTime,~,XP,TP,~,~,~,W] = Cheb_prediction(PredictTime,recPointsT,Xtrain,Dtrain,hT,L,pM,CC,bet,Xpred,Tpred,plt);


pM2 = linspace(0.58,0.75,150);
XpredMaxArray=[];
XsignalMaxArray=[];
for i1=1:length(pM2)
    P.b =pM2(i1);
    
    [X2,~]=Generate_signal(hT,NT,x0,P);
    x0 = X2(:,end);
    [X3,~]=Generate_signal(hT,NT,x0,P);
    x0 = X3(:,end);



    if i1 ==1
        LL = length(X2(1,L-2*Ntau:end) );
        Xpred = zeros(dim+dim_c,LL);
        X33 = X3;
        X3 = zeros(dim,NT);
        for i2 = 1:dim
            tmpX = X2(varbl(i2),L-2*Ntau:end);
            Xpred(i2,:) = resc(i2)*(tmpX-Xmn(i2))-1;
            tmpX3 = X33(varbl(i2),:);
            X3(i2,:) = resc(i2)*(tmpX3-Xmn(i2))-1;
        end
    else
        X33 = X3;
        X3 = zeros(dim,NT);
        for i2 = 1:dim
            tmpX3 = X33(varbl(i2),:);
            X3(i2,:) = resc(i2)*(tmpX3-Xmn(i2))-1;
        end
    end

    LL = length(Xpred(1,:) );
    Xpred(dim+1,:) = ones(1,LL)*P.b;



    XP=integrate_reconstructed_system(W,Xpred,hT,NT,CC,Ntau,deg);


    Xpred= XP(:,end-(k*dim-1)*Ntau-1:end);
     
    XP=integrate_reconstructed_system(W,Xpred,hT,NT,CC,Ntau,deg);

    jX3 = islocalmax(X3(1,:));
    jXP = islocalmax(XP(1,:));

    sXP ='r.';
    sX3 ='b.';

    %if solution is stationary
    if isempty(find(jXP,1))
        jXP = length(XP(1,:));
        sXP = 'r*';
    end

    if isempty(find(jX3,1))
        jX3 = length(X3(1,:));
        sX3 = 'b*';
    end

    jmX3 = islocalmin(X3(1,:));
    jmXP = islocalmin(XP(1,:));
    
    tmpX3 = 1:length(X3(1,:));
    tmpXP = 1:length(XP(1,:));

    XPP=XP;
    for i2 =1:dim
        tmpXP = XP(i2,:);
        XPP(i2,:) = (tmpXP+1)/resc(i2)+Xmn(i2);
    end

%     figure(79)
%     subplot(311)
%     ss=scatter(X3(1,jX3)*0+P.b,X33(varbl(1),jX3),sX3,'SizeData',5); hold on
%     ss.MarkerEdgeAlpha = 0.5;
%     subplot(312)
%     ss=scatter(XP(1,jXP)*0+P.b,XPP(1,jXP),sXP,'SizeData',5); hold on
%     ss.MarkerEdgeAlpha = 0.5;
%     subplot(313)
%     ss1=scatter(X3(1,jX3)*0+P.b,X3(1,jX3),sX3,'SizeData',5); hold on
%     ss2=scatter(XP(1,jXP)*0+P.b,XP(1,jXP),sXP,'SizeData',5); hold on
%     ss1.MarkerEdgeAlpha = 0.5;
%     ss2.MarkerEdgeAlpha = 0.5;

    XpredMaxArray = [XpredMaxArray;XPP(1,jXP)'*0+P.b,XPP(1,jXP)' ];
    XsignalMaxArray = [XsignalMaxArray;X33(varbl(1),jX3)'*0+P.b,X33(varbl(1),jX3)' ];
    
    
end

% figure(79)
% subplot(311)
% ylabel('actual')
% subplot(312)
% ylabel('prediction')

save('results2.mat','XsignalMaxArray','XpredMaxArray',"pM2","pM")
