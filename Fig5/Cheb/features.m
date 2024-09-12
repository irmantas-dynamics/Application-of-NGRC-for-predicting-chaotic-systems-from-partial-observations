function R = features(Cheb,X,polyn_deg,dtrvec,dtot,szL)

    dpol = length(polyn_deg);
    R=zeros(dtot,szL);
    idx1 = 1;
    idx2 = 0;
    for i1 =1:dpol
        idx2= idx2+dtrvec(i1);
        tmp_train=Cheb(X,polyn_deg(i1),dtrvec(i1));
        R(idx1:idx2,:)= tmp_train;
        idx1=idx2+1;
    end

end