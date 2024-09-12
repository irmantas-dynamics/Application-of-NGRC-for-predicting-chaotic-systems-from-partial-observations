function R=prepare_feature_vector(X0,L,Ntau,k,dim,dim_c,dpol,polyn_deg)  
    szL = ceil(L-(k-1)*Ntau);
    
    X=zeros(k*dim+dim_c,szL);    
    count = 0;
    for i1 = 1:k
        for i2 = 1:dim
            count= count+1;
            buf = (k-i1)*Ntau;
            buf2 = (i1-1)*Ntau;
            X(count,:)= X0(i2,buf+1:L-buf2);
        end
    end

    for i1= 1:dim_c
        count= count+1;
        X(count,:) = X0(i1+dim,(k-1)*Ntau+1:L);
    end

 

    

    % vector dtrvec(i1) save sizes of vector combined from all possible
    % polynomials with degree polyn_deg(i1)
    dtrvec = ones(dpol,1);
    for i1 =1:dpol
        for i2 =1:polyn_deg(i1)
            dtrvec(i1) = dtrvec(i1)*((k*dim+dim_c)+i2-1)/i2;
        end
    end
    dtot = sum(dtrvec);

    %Chebyshev n-order, for k*dim+dim_c variables
    fnname = sprintf('Cheb_Rn_%d',k*dim+dim_c);
    if ~exist(strcat(fnname,'.m'),'file')
        generate_Polyn(k*dim+dim_c);
    end

    Cheb = str2func(fnname);

    R =features(Cheb,X,polyn_deg,dtrvec,dtot,szL);
end