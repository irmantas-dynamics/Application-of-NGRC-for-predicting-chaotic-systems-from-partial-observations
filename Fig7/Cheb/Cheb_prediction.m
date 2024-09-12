function [Error,ErrorT,ErrorL,countNaN,XP,TP,R,RR,D,W]=Cheb_prediction(PredictTime,X0,hT,CC,bet,Xpred,Tpred,plt)
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

    L=length(X0)-1; % Tasku skaicius apmokymui
        
    k =  CC.k; % number of embeded coordinates
    
    % polynomial degrees that would be used for regression
    polyn_deg = CC.polyn_deg;
    dpol = length(polyn_deg);

%     Xmn = min(X0);
%     Xmx = max(X0);
%     a = 2/(Xmx-Xmn);
%     X0 = a*(X0-Xmn)-1;
    M=max(max(X0),-min(X0));
    X0 =X0/M;

    szL = ceil((L-(k-1)*Ntau)/nT);
    
    X=zeros(k,szL);    
    for i1 = 1:k
        buf = (k-i1)*Ntau;
        buf2 = (i1-1)*Ntau;
        X(i1,:)= X0(buf+1:nT:L-buf2);
    end
    

    % vector dtrvec(i1) save sizes of vector combined from all possible
    % polynomials with degree polyn_deg(i1)
    dtrvec = ones(dpol,1);
    for i1 =1:dpol
        for i2 =1:polyn_deg(i1)
            dtrvec(i1) = dtrvec(i1)*(k+i2-1)/i2;
        end
    end
    dtot = sum(dtrvec);

    %Monomials n-order, for k variables
    fnname = sprintf('Cheb_Rn_%d',k);
    if ~exist(strcat(fnname,'.m'),'file')
        generate_Polyn(k);
    end

    Cheb = str2func(fnname);

    R =features(Cheb,X,polyn_deg,dtrvec,dtot,szL);

    D1 = (X0(2:L+1)-X0(1:L))/h;
    D=D1((k-1)*Ntau+1:nT:L);

    d=size(R,1);
    RR=R*R.'+bet*eye(d);
    B=R*D.';
    W = (RR\B).';
%     W = (pinv(RR)*B).';
    ErrorL = sqrt(mean((D-W*R).^2))/std(D);
    if plt ==1
        fprintf('====================================\n\n')
        fprintf('\t\t\t PREDICTION TIME: %f\n\n',PredictTime/CC.lyapTime)
        fprintf('Learning error: %f\n',ErrorL)
        fprintf('Number of unknown variables: %d\n',dtot)
        fprintf('Number of feature vectors: %d\n',szL)

%         sqrt(mean((D-W*R).^2))
    end



    % % Prognoze po apmokymo
%     Tp=(k-1)*Ntau*h+ round(PredictTime/h)*h; %predikcijos laikas

%     Y1 =  a*(Xpred-Xmn)-1;
    Y1 = Xpred /M;
   
    Z11 = Y1;
 
    Np=round(PredictTime/h)+(k-1)*Ntau+1; % Predikcijos masyvo ilgis
    NAllP = round(Tpred(end)/h)-Np+2;


    TP=Tpred(1:Np);
    XP=NaN(1,Np);

    XPiter = zeros(NAllP,1);

    countNaN = 0;


    for i2 =1:NAllP
        % Oilerio integravimas rekonstruotos DDE
        
        XP(1:(k-1)*Ntau+1)= Z11(i2:(k-1)*Ntau+i2);
        x=zeros(k,1);    
        for j=(k-1)*Ntau+1:Np-1
            for i1 = 1:k
                buf = (i1-1)*Ntau;
                x(i1)= XP(j-buf);
            end
            r =features(Cheb,x,polyn_deg,dtrvec,dtot,1);
            XP(j+1)=x(1)+h*(W*r);
            if abs(XP(j+1)) >2
                 if (Z11(Np+i2:end,1)-1)<-1
                    XP(Np) = 1;
                 else
                    XP(Np) =-1;
                 end

                 countNaN = countNaN+1;
                break;
            end
        end
        XPiter(i2)=XP(Np);
    end

    Error  = sqrt(mean( (Z11(Np:end)-XPiter').^2))/std(Z11(Np:end));
    ErrorT = Z11(Np:end)-XPiter';

   
    if plt == 1
        figure
        subplot(211)
        plot(Tpred(Np:end),Z11(Np:end)), hold on
        plot(Tpred(Np:end),XPiter,'--')
        legend({'True signal','prediction'})
        title(sprintf('Pred time(lyap): %f'),PredictTime/CC.lyapTime)
        subplot(212)
        plot(Tpred(Np:end),Z11(Np:end)-XPiter')
        fprintf('STD of true signal: %f\n',std(Z11(Np:end)));
        fprintf('Prediction Error: %f\n',Error);
        fprintf('====================================\n\n')
    end

end
 
