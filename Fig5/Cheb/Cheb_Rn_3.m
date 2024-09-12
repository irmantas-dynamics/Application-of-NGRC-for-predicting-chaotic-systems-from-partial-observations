function r=Cheb_Rn_3(X,n,N)

	[k, L]=size(X);
	x1=zeros(n+1,L);
	x2=zeros(n+1,L);
	x3=zeros(n+1,L);

	x1(1,:)=ones(1,L);
	x2(1,:)=ones(1,L);
	x3(1,:)=ones(1,L);

	x1(2,:) = X(1,:);
	x2(2,:) = X(2,:);
	x3(2,:) = X(3,:);

	for i=2:n
		x1(i+1,:)=2*x1(i,:).*X(1,:)-x1(i-1,:);
		x2(i+1,:)=2*x2(i,:).*X(2,:)-x2(i-1,:);
		x3(i+1,:)=2*x3(i,:).*X(3,:)-x3(i-1,:);
	end

	r=zeros(N,L);
	m=0;
	for k1=0:n
		for k2=0:n-k1
			k3= n-k1-k2;
			m=m+1;
			r(m,:)= x1(k1+1,:).*x2(k2+1,:).*x3(k3+1,:);
		end
	end
end
