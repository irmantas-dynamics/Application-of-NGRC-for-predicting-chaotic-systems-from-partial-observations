function run_program_function(deg,betdeg)

    recPoints = [0.1; 0.5; 0.75;1;1.5]; % Errors recording points at time scale
    hTM= [0.001; 0.005; 0.05; 0.15];
    % hTM = [0.03; 0.05];
    Ttrain = 200; %training time
    
    
    load("init.mat","x0M");
    szx0=size(x0M,1);
%     szx0 = 20; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % which initial condition from x0M array to use
    % iix = 10;
    
    
    
    i2 = betdeg; % i2    - abs() of beta degree
    bet =  10^(-i2);
    
    CC.lyapTime = 1.104;
    
    plt = 0; % turn on (1)/off(0) plotting
    
    
    % Total prediction time
    PredictTime = 10*CC.lyapTime;   
    
    % deg = 1:2:15;
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
    
    
    trainErrM = NaN(length(hTM));
    predErrM  = NaN(length(hTM),length(recPoints));
    nanCountM = zeros(length(hTM),length(recPoints));
    MeanTime = zeros(length(hTM),1);
    STDTime = zeros(length(hTM),1);
    quantileTime = zeros(length(hTM),3);
    
    
    
    histSteps = 0:0.1:10;
    histM1 = zeros(length(hTM),length(histSteps)-1);
    
    recPointsT  =recPoints*CC.lyapTime ;
    
    tic
    for ihT=1:length(hTM)
    %     disp ('==========================================')
        hT =hTM(ihT);
        NT = ceil((Ttrain+PredictTime)/hT)+1;
        L = ceil(Ttrain/hT)+1;
    
        SpredErr = zeros(szx0,length(recPoints));
        StrainErr = 0;
        SnanCount = zeros(1,length(recPoints));
        Terr = zeros(1,szx0);
        
    %     pip = 10;
        parfor ix0= 1:szx0
            x0=x0M(ix0,:);
            [X1,~]=Generate_signal(hT,NT,x0,P);
             T=linspace(0,hT*(NT-1),NT);
    
    %         figure(1)
    %         plot(T,X1(1,:),'.'), hold on
    
            X0 = X1(varbl,1:L);
            Ntau=round(tau/hT);
            Xpred = X1(varbl,L-2*Ntau:end);
            Tpred = T(L-2*Ntau:end)-T(L-2*Ntau);
    
            [ErrorL,Terr(1,ix0),ErrorFixTime,~,XP,~,~,~,~,~] = Cheb_prediction(PredictTime,recPointsT,X0,hT,CC,bet,Xpred,Tpred,plt);
            StrainErr =StrainErr + ErrorL;
            nnan = isnan(XP)';
            SnanCount=SnanCount+nnan;
    
    
        end
            [histM1(ihT,:),edges]=histcounts(Terr(1,:),histSteps);
            trainErrM(ihT) = StrainErr/szx0;
            MeanTime(ihT) = mean(Terr(1,:));
            STDTime(ihT) = std(Terr(1,:));
            quantileTime(ihT,:) = quantile(Terr(1,:),[0.25,0.5,0.75]);
    end
    
    fpav = sprintf('Cheb%d_bet%d.mat',max(deg),i2);
    save(fpav)

end
