function run_program_function(deg,betdeg)
    
    
    % deg = 1:2:21;
    i2 = betdeg; % i2    - abs() of beta degree
    bet =  10^(-i2);
    
    
    hT =  0.05;
    Ttrain = 1000; %training time
    
    
    load("init.mat","x0M");
    szx0=size(x0M,1);
    szx0 = 24;
    
    
    % % which initial condition from x0M array to use
    % iix = 10;
    
    CC.lyapTime = 1.104;
    
    plt = 0; % turn on (1)/off(0) plotting
    
    
    % Prediction time
    % PredictTime = 0.7*CC.lyapTime;  
    PredictTimeM = [0.1;0.3;0.5;0.7;0.9]*CC.lyapTime;  
    % PredictTimeM = [1.]*CC.lyapTime;  
    
    %total prediction time 
    Tpred = 400;
    
    
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
    
    
    trainErrM = zeros(szx0,length(PredictTimeM));           
    ErrorM = zeros(szx0,length(PredictTimeM));
    countNaNM = zeros(szx0,length(PredictTimeM));
    
    
    
    %     disp ('==========================================')
        NT = ceil((Ttrain+Tpred)/hT)+1;
        L = ceil(Ttrain/hT)+1;
       
    %     pip = 1;
        parfor ix0= 1:szx0
    %     for ix0 = pip:pip
            x0=x0M(ix0,:);
            [X1,~]=Generate_signal(hT,NT,x0,P);
            T=linspace(0,hT*(NT-1),NT);
    
            X0 = X1(varbl,1:L);
            Ntau=round(tau/hT);
            Xpred = X1(varbl,L-2*Ntau:end);
            Tpred = T(L-2*Ntau:end)-T(L-2*Ntau);
    
            trainErrtmp = zeros(1,length(PredictTimeM));
            Errortmp = zeros(1,length(PredictTimeM));
            countNaNtmp = zeros(1,length(PredictTimeM));
    
            for iPt =1:length(PredictTimeM)
                PredictTime = PredictTimeM(iPt); 
        
                [Errortmp(iPt),~,trainErrtmp(iPt),countNaNtmp(iPt),XP,~,~,~,~,~] = Monom_prediction(PredictTime,X0,hT,CC,bet,Xpred,Tpred,plt);
    
            end
            trainErrM(ix0,:)= trainErrtmp;
            ErrorM(ix0,:) = Errortmp;
            countNaNM(ix0,:) = countNaNtmp;
    
    
        end
    
    
    fpav = sprintf('Monom%d_bet%d.mat',max(deg),i2);
    save(fpav)
    
end

