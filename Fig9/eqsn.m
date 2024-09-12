function yp = eqsn(~,y,P)
    % El_Nino_Southern sistema
    B=P.B;
    Delx=P.Delx;
    C=P.C;
    u0=P.u0;
    Tm=P.Tm;
    A=P.A;
    T0=P.T0;

    u=y(1);
    Tw=y(2);
    Te=y(3);
    
    yp=[B*(Te-Tw)/(2*Delx)-C*(u-u0);
        u*(Tm-Te)/(2*Delx)-A*(Tw-T0);
        u*(Tw-Tm)/(2*Delx)-A*(Te-T0)];
end