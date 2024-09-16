function run_program_function(deg,betdeg)
% deg - degree of polynome
% betdeg - abs() of beta degree
    
    %reguliarization parameter
    bet =  10^(-betdeg);
    
    %discretization step
    hT =  0.1;
    %training time
    Ttrain = 2000; %training time

    %loading initial conditions
    load("init.mat","x0M");
    % we will use only the first 24 conditions
    szx0 = 24;

    CC.lyapTime = 14.08;

    plt = 0; % turn on (1)/off(0) intermediate plotting


    % Continuous prediction time
    PredictTimeM = [2;4;6;8;10]*CC.lyapTime;  

    %total prediction time 
    Tpred = 40*CC.lyapTime;


    %number of embedded dimension
    k= 3;

    % delay of embeded dimensions
    tau = 1.5;

    CC.tau= tau;
    CC.k = k;
    CC.polyn_deg = deg;
    %which variable will be predicted
    varbl = 2;

    %Model (Rossler) parameters
    a= 0.2;
    b= 0.2;
    c= 5.7;

    P.a=a;
    P.b=b;
    P.c=c;

    %array of the training error
    trainErrM = zeros(szx0,length(PredictTimeM));           
    %array of the prediction error
    ErrorM = zeros(szx0,length(PredictTimeM));
    %array of how many trajectories from the given inital conditions diverge to infinity
    countNaNM = zeros(szx0,length(PredictTimeM));


    NT = ceil((Ttrain+Tpred)/hT)+1;
    L = ceil(Ttrain/hT)+1;
   
    parfor ix0= 1:szx0 %use this if you have parralel computing toolbox, comment otherwise
%     for ix0 = 1:szx0 %uncoment if you have not parralel computing toolbox
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
           [Errortmp(iPt),~,trainErrtmp(iPt),countNaNtmp(iPt),XP,~,~,~,~,~] = Cheb_prediction(PredictTime,X0,hT,CC,bet,Xpred,Tpred,plt);

        end
        trainErrM(ix0,:)= trainErrtmp;
        ErrorM(ix0,:) = Errortmp;
        countNaNM(ix0,:) = countNaNtmp;
    end

    fpav = sprintf('Cheb%d_bet%d.mat',max(deg),betdeg);
    save(fpav)
end

