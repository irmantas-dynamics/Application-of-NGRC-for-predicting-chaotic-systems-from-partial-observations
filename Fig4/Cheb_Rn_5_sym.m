function r=Cheb_Rn_5_sym(X,n,N)

	[k, L]=size(X);
	x1=sym('a',[n+1,1]);
	x2=sym('a',[n+1,1]);
	x3=sym('a',[n+1,1]);
	x4=sym('a',[n+1,1]);
	x5=sym('a',[n+1,1]);

	x1(1,:)=ones(1,L);
	x2(1,:)=ones(1,L);
	x3(1,:)=ones(1,L);
	x4(1,:)=ones(1,L);
	x5(1,:)=ones(1,L);

	x1(2,:) = X(1,:);
	x2(2,:) = X(2,:);
	x3(2,:) = X(3,:);
	x4(2,:) = X(4,:);
	x5(2,:) = X(5,:);

	for i=2:n
		x1(i+1,:)=x1(i,:).*X(1,:);
		x2(i+1,:)=x2(i,:).*X(2,:);
		x3(i+1,:)=x3(i,:).*X(3,:);
		x4(i+1,:)=x4(i,:).*X(4,:);
		x5(i+1,:)=x5(i,:).*X(5,:);
	end

	r=sym('a',[N,1]);
	m=0;
	for k1=0:n
		for k2=0:n-k1
			for k3=0:n-k1-k2
				for k4=0:n-k1-k2-k3
					k5= n-k1-k2-k3-k4;
					m=m+1;
					r(m,:)= x1(k1+1,:).*x2(k2+1,:).*x3(k3+1,:).*x4(k4+1,:).*x5(k5+1,:);
				end
			end
		end
	end
end
