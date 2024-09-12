function XP=integrate_reconstructed_system(W,x0,h,Np,CC,Ntau,polyn_deg)

    k= CC.k;
    dim = CC.dim;
    dim_c =CC.dim_c;
    dpol = length(polyn_deg);
    

    %Monomials n-order, for k variables
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


    XP=NaN(dim+dim_c,Np);

    % Oilerio integravimas rekonstruotos DDE
    XP(:,1:(k*dim-1)*Ntau+1)= x0(:,1:(k*dim-1)*Ntau+1);
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

end