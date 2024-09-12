function  [Y1,DY1]=Generate_signal(hT,NT,y0,P)
 
    opts = ddeset('RelTol',1e-7,'AbsTol',1e-10);
    
    T=linspace(0,hT*(NT-1),NT);
    %[~,Y1] = ode45(@sistema,T,y00,opts);
    [~,Y1] = ode45(@(t,y)eqsn(t,y,P), T,y0,opts);
    
    Y1= Y1.';   
    DY1 = [];

end
%---------------------------------------------------------------------