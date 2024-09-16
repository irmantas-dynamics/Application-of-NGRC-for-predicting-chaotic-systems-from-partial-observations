function generate_Polyn(k)
    
    pav = sprintf('Cheb_Rn_%d.m',k);
    fid = fopen(pav,'w');
    
    pav = sprintf('function r=Cheb_Rn_%d(X,n,N)\n\n',k);
    fprintf(fid,pav);

    fprintf(fid,"\t[k, L]=size(X);\n");
%     fprintf(fid,"\tN=round((n+1)*(n+2)/2);\n");

    for k1=1:k
        pav = sprintf("\tx%d=zeros(n+1,L);\n",k1);
        fprintf(fid,pav);
    end

    fprintf(fid,"\n");

    for k1=1:k
        pav = sprintf("\tx%d(1,:)=ones(1,L);\n",k1);
        fprintf(fid,pav);
    end

     fprintf(fid,"\n");




    for k1=1:k
        pav = sprintf("\tx%d(2,:) = X(%d,:);\n",k1,k1);
        fprintf(fid,pav);
    end

    fprintf(fid,"\n");
    fprintf(fid,"\tfor i=2:n\n");

    for k1=1:k
        pav = sprintf("\t\tx%d(i+1,:)=2*x%d(i,:).*X(%d,:)-x%d(i-1,:);\n",k1,k1,k1,k1);
        fprintf(fid,pav);
    end
    fprintf(fid,"\tend\n\n");

    fprintf(fid,"\tr=zeros(N,L);\n");
    fprintf(fid,"\tm=0;\n");

    
    for i1 =1:k-1
        pav =sprintf('for k%d=0:n',i1);
        pav =strcat(repmat('\t', 1, i1),pav);
        for i2 = 1:i1-1
            pav2 =sprintf('-k%d',i2);
            pav =strcat(pav,pav2);
        end
        pav =strcat(pav,'\n');
        fprintf(fid,pav);
    end

    pav = sprintf('k%d= n',k);
    pav =strcat(repmat('\t', 1, k),pav);
    for i2 = 1:k-1
       pav2 =sprintf('-k%d',i2);
       pav =strcat(pav,pav2);
    end
    pav =strcat(pav,';\n');
    fprintf(fid,pav);

    pav = 'm=m+1;\n';
    pav =strcat(repmat('\t', 1, k),pav);
    fprintf(fid,pav);

    pav = 'r(m,:)= x1(k1+1,:)';
    for i2 = 2:k
       pav2 =sprintf('.*x%d(k%d+1,:)',i2,i2);
       pav =strcat(pav,pav2);
    end
    pav =strcat(pav,';\n');
    pav =strcat(repmat('\t', 1, k),pav);
    fprintf(fid,pav);

    for i1 = k-1:-1:1
        fprintf(fid,strcat(repmat('\t', 1, i1),'end\n'));

    end

% 
%     fprintf(fid,"\tfor i=0:n\n");
%     fprintf(fid,"\t\tfor j=0:n-i\n");
%     fprintf(fid,"\t\t\tk=n-i-j;\n");
%     fprintf(fid,"\t\t\tm=m+1;\n");
%     fprintf(fid,"\t\t\tr(m,:)=x(i+1,:,1).*x(j+1,:,2).*x(k+1,:,3);\n");
%     fprintf(fid,"\t\tend\n");
%     fprintf(fid,"\tend\n");
    fprintf(fid,"end\n");

    fclose(fid);

end