clear all
A=[10,7,8,7;7,5,6,5;8,6,10,9;7,5,9,10];%这只是个例子，你可以任意改变它
fprintf('matlab函数计算的特征值为： %.8f,%.8f,%.8f,%.8f\n',eig(A))
%hess(A)
%Householder transformations to Hessenberg matrix
n=length(A(:,1));
for i=2:n
    sigma=sign(A(i,i-1))*sqrt(sum(A(i:n,i-1).^2));
    beta=sigma*(A(i,i-1)+sigma);
    A(i,i-1)=A(i,i-1)+sigma;
    R=eye(n-i+1)-1.0/beta*A(i:n,i-1)*A(i:n,i-1)';
    A(i,i-1)=-sigma;
    A(i+1:n,i-1)=0;
    A(i-1,i:n)=A(i:n,i-1);%AT=A
    %A(1:i-1,i:n)=A(1:i-1,i:n)*R;%AT\=A
    A(i:n,i:n)=R*A(i:n,i:n)*R;
end
%upper Hessenberg Matrix QR
%Givens transformation
H=A;
m=n;
c=zeros(n,1);
d=zeros(n,1);
while (m>2)
    s=H(m,m);
    H(1,1)=H(1,1)-s;
    for k=1:m-1
        H(k+1,k+1)= H(k+1,k+1)-s;
        cs=sqrt(H(k,k)^2+H(k+1,k)^2);
        c(k)=H(k,k)/cs;d(k)=H(k+1,k)/cs;
        for i=k:m
           tmp=[c(k),d(k);-d(k),c(k)]*[H(k,i),H(k+1,i)]';
           H(k,i)=tmp(1);H(k+1,i)=tmp(2);
        end
    end
    for k=1:m-1
        for i=1:k+1
            tmp=[H(i,k),H(i,k+1)]*[c(k),-d(k);d(k),c(k)];
             H(i,k)=tmp(1);H(i,k+1)=tmp(2);
        end
        H(k,k)=H(k,k)+s;
    end
    H(m,m)=H(m,m)+s;
    tol=abs(H(m,m-1));
    while(tol<0.00000000001 & m>2)
        m=m-1;
        tol=abs(H(m,m-1));
    end
end
for i=1:n
    fprintf('A的第%d个特征值为：%0.8f\n',i,H(i,i))
end
