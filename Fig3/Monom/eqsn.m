function dy = eqsn(~,y,P)
%rossler system
    a= P.a;
    b= P.b;
    c= P.c;

    dy=zeros(3,1);
    dy(1) = -y(2)-y(3);
    dy(2) = y(1)+a*y(2);
    dy(3) = b+y(3)*(y(1)-c);
end