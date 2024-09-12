function [ErrorL,Terr,ErrorFixTime,signal,XP,TP,R,RR,D,W]=Cheb_prediction(PredictTime,recPoints,Xtrain,Dtrain,hT,L,pM,CC,bet,Xpred,Tpred,plt)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Cheb_prediction() from signal X0 makes prediction by using NVAR with
% Chebyshev polynomes. Function also compares its predictions with signal
% from orginal system Xpred.
%
% Input:
%   PredictTime - prediction time, scalar
%   recPoints - time points to record error,array
%   X0 -  Signal for learning process, array
%   hT - integration step
%   CC - learning parameters, structure
%       CC.k - embeding dimension size
%       CC.polyn_deg - degrees of polynomes, that will be used in featrue
%                      vector, vector. Eg. [1,2,5,7] means that
%                      x,x^2,x^5,x^7 will be used in approximation
%       CC.tau - minimal delay of embeded dimension, scalar
%   bet - reguliarization parameter, scalar
%   Xpred - signal that shoul be predicted, vector
%   plt - ploting option, 1 - plots figures; 0 -does not plot, scalar
%
% Output: 
%   ErrorL - learning error, scalar
%   Terr - Time when prediction deviates from true solution. (Threshold hardcoded)
%   ErrorFixTime - difference between true and predicted signals at given
%                  (recPoints) time moments, vector
%   signal - returns Xpred, vector
%   XP -  predicted signal, includes initial conditions for delayed
%         values, i.e -tau...0..., vector
%   TP - predicted signal time, negative TP values used as initial
%        conditions, for prediction algorithm, vector
%   R - feature vector
%   RR  - R*R.'+bet, matrix
%   D - derivative of X0, vector 
%   W - learned matrix
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    tau = CC.tau;
    h=hT; % laiko zingsnis
    Ntau=round(tau/h); % turi buti sveikas skaicius
    nT= 1;% kas kelinta zingsni imti isvestines fitinimui   

%     L=size(X0,2)-1; % Tasku skaicius apmokymui
        
    k =  CC.k; % number of embeded coordinates
    dim = CC.dim;
    %dim_c number of constant inputs
    dim_c = CC.dim_c;
    % Cmn(i) and Cmx(i)  - range of i-th constant input (e.g. for bifurcation diagram we record the dynamics of variable at several constant parameter values)
    % length(Cmn or Cmx) = dim_c
    Cmn = CC.Cmn;
    Cmx = CC.Cmx;
    
    % polynomial degrees that would be used for regression
    polyn_deg = CC.polyn_deg;
    dpol = length(polyn_deg);


    R=[];
    D=[];
    for i1 = 1:length(pM)
        X0 = Xtrain(:,1+(i1-1)*L:i1*L);
        D0 = Dtrain(:,1+(i1-1)*L:i1*L);
        R0 = prepare_feature_vector(X0,L,Ntau,k,dim,dim_c,dpol,polyn_deg);
        R=[R R0];
        D1=D0(:,(k-1)*Ntau+1:L);
        D= [D D1];
    end



    d=size(R,1);
    RR=R*R.'+bet*eye(d);
    B=R*D.';
    W = (RR\B).';
%     W = (pinv(RR)*B).';
    ErrorL = sqrt(mean((D-W*R).^2))/std(D);
    if plt ==1
        disp('Learning error')
        disp(ErrorL)
%         sqrt(mean((D-W*R).^2))
    end


    %Cheb n-order, for k variables
    fnname = sprintf('Cheb_Rn_%d',k*dim+dim_c);
    if ~exist(strcat(fnname,'.m'),'file')
        generate_Polyn(k*dim+dim_c);
    end

    Cheb = str2func(fnname);


    % vector dtrvec(i1) save sizes of vector combined from all possible
    % polynomials with degree polyn_deg(i1)
    dtrvec = ones(dpol,1);
    for i1 =1:dpol
        for i2 =1:polyn_deg(i1)
            dtrvec(i1) = dtrvec(i1)*((k*dim+dim_c)+i2-1)/i2;
        end
    end
    dtot = sum(dtrvec);


    % % Prognoze po apmokymo
    Tp=(k-1)*Ntau*h+ round(PredictTime/h)*h; %predikcijos laikas


    Z11 = Xpred;

    TP=Tpred;
    Np=length(TP); % Predikcijos masyvo ilgis
    XP=NaN(dim+dim_c,Np);


    % Oilerio integravimas rekonstruotos DDE
    XP(:,1:(k*dim-1)*Ntau+1)= Z11(:,1:(k*dim-1)*Ntau+1);
    x=zeros(dim*k+dim_c,1);  
    for j=(k*dim-1)*Ntau+1:Np-1

        count =0;
        for i1 = 1:k
            for i2=1:dim
                count = count+1;
                buf = (i1-1)*Ntau;
                x(count)= XP(i2,j-buf);
            end
        end

        for i1= 1:dim_c
            count= count+1;
            x(count) = XP(i1+dim,j);
        end

        r =features(Cheb,x,polyn_deg,dtrvec,dtot,1);
        XP(1:dim,j+1)=x(1:dim,1)+h*(W*r);
        XP( (dim+1):(dim_c+dim),j+1)= XP((dim+1):(dim_c+dim),j);
    end

    idx = (k-1)*Ntau+1;
    TP=TP-TP(idx);


    T1= TP;
    lnT = length(T1);

    skirt = abs(Z11(1:Np)-XP)/2;

    idxRec = round(recPoints/h);
    ErrorFixTime = zeros(length(idxRec),1);

    if plt == 1
        figure
        subplot(311)
        plot(TP/CC.lyapTime,XP(1,:)), hold on
        plot(TP/CC.lyapTime,Z11(1,1:lnT),'--')
        plot(TP(idx)/CC.lyapTime,XP(1,idx),'r*')
%         ylim([-1.1 1.1])
        xlim([0 max(TP/CC.lyapTime)])
        legend({'Predict','real'})

        subplot(312)
        plot(TP/CC.lyapTime,XP(2,:)), hold on
        plot(TP/CC.lyapTime,Z11(2,1:lnT),'--')
        plot(TP(idx)/CC.lyapTime,XP(2,idx),'r*')
%         ylim([-1.1 1.1])
        xlim([0 max(TP/CC.lyapTime)])


%         subplot(313)
%         plot(TP/CC.lyapTime,XP(3,:)), hold on
%         plot(TP/CC.lyapTime,Z11(3,1:lnT),'--')
%         plot(TP(idx)/CC.lyapTime,XP(3,idx),'r*')
% %         ylim([-1.1 1.1])
%         xlim([0 max(TP/CC.lyapTime)])

        figure
        imagesc(W)

        fnname = sprintf('Cheb_Rn_%d_sym',k*dim+dim_c);
        if ~exist(strcat(fnname,'.m'),'file')
            generate_Polyn_sym(k*dim+dim_c);
        end
    
%         fun = str2func(fnname);

%         syms x1 x2 x3 p
%         Olin_i = [x1;x2;x3;p];
% %         Olin_i = [x_1;y_1;z_1;  p];
% 
%         rsym =features_sym(fun,Olin_i,polyn_deg,dtrvec,dtot,1);
% 
%         
%         xticks(1:length(rsym))
% %         xticklabels({Otot_i})
% 
%         xtik = {};
%         for i1 = 1:length(rsym)
%             tmp = sprintf('%s',rsym(i1));
%             tmp = erase(tmp,'*');
%             xtik{i1} = tmp;
%         end
% 
% 
%         xticklabels(xtik)

%         subplot(212)
%         plot(TP/CC.lyapTime,skirt), hold on
%         plot(TP/CC.lyapTime,TP/CC.lyapTime*0+0.01,'r-')
%         ylim([-0.1 2.1])
%         xlim([0 max(TP/CC.lyapTime)])
%         ylabel('Error')
%         xlabel('\Lambda_{max} t')

    end

%     for i1=1:length(idxRec)
%         if plt ==1 
%             subplot(211)
%             plot(TP(idx+idxRec(i1))/CC.lyapTime , XP(idx+idxRec(i1)),'ro')
%             plot(T1(idx+idxRec(i1))/CC.lyapTime , Z11(idx+idxRec(i1)),'bx')
%         end
%         idx2 = idx+idxRec(i1);
%         ErrorFixTime(i1) = abs( XP(idx2)- Z11(idx2));
%     end


%     save('duom2.mat')

    ind1 = find(skirt>0.01);
    if isempty(ind1)
        ind1 = lnT;
    else
        ind1 =ind1(1);
    end  

    Terr = T1(ind1);

    signal = Z11(Np:end);

end
 

