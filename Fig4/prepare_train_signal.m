function [Xtrain, Dtrain,a,Xmn] = prepare_train_signal(hT,L,x0,P,pM,varbl)

    Np = length(pM);
    dim= length(varbl);
    Xtrain = zeros(dim+1,Np*L);
    Dtrain = zeros(dim,Np*L);
    
    x00=x0;
    for i1=1:Np
        P.b = pM(i1);

        [X1,~]=Generate_signal(hT,L,x00,P);
        x00=X1(:,end);

        [X1,~]=Generate_signal(hT,L+1,x00,P);
        x00=X1(:,end);
    
        D1 = (X1(varbl,2:L+1)-X1(varbl,1:L))/hT;
        Xtrain(:,1+(i1-1)*L:i1*L)=[X1(varbl,1:L); ones(1,L)*pM(i1)];
        Dtrain(:,1+(i1-1)*L:i1*L)=D1;
    end

    Xmn = zeros(dim,1);
    a = zeros(dim,1);
    for i1 = 1 :dim
        Xmn(i1) = min(Xtrain(i1,:));
        Xmx = max(Xtrain(i1,:));
        a(i1) = 2/(Xmx-Xmn(i1));
        Xtrain(i1,:) = a(i1)*(Xtrain(i1,:)-Xmn(i1))-1;
        Dtrain(i1,:) = a(i1)*Dtrain(i1,:);
    end


end