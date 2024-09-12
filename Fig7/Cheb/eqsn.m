function yp = eqsn(~,y,P)
%lorenz system
sig=P.sig; 
r=P.r;
b=P.b;
yp=[-sig*(y(1)-y(2));
     r*y(1)-y(2)-y(1)*y(3);
     y(1)*y(2)-b*y(3)];
end
