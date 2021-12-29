clear all
A=[6,2,1;2,3,1;1,1,1];
s=6;
n=length(A(:,1));
u0=ones(n,1);%初始向量
v0=zeros(n,1);
%=====LU=求解Uv0=u0===
[B,p]=LUp(A-s*eye(n));
v0(n)=u0(n)/B(n,n);
for i=n-1:-1:1
    tmp=0;
    for j=i+1:n
        tmp=tmp+B(i,j)*v0(j);
    end
    v0(i)=(u0(i)-tmp)/B(i,i);
end
tol=1;
mu1=max(abs(v0));
u0=v0/mu1;
total =1;
while(tol >0.0000001)
    a=mu1;
    %====Pu0 换行=======
    for i=1:n
        if(p(i)==i)
            continue;
        end
            tmp1=u0(i);
            u0(i)=u0(p(i));
            u0(p(i))=tmp1;
    end
    %===Ly=Pu0============
    for i=2:n
        tmp=0;
        for k=1:i-1
            tmp=tmp+B(i,k)*u0(k);
        end
        u0(i)=u0(i)-tmp;
    end
    %====Ux=y===========
    v0(n)=u0(n)/B(n,n);
    for i=n-1:-1:1
        tmp=0;
        for j=i+1:n
            tmp=tmp+B(i,j)*v0(j);
        end
        v0(i)=(u0(i)-tmp)/B(i,i);
    end
   mu1=max(abs(v0));
   u0=v0/mu1;
   total =total +1;
   tol=abs(mu1-a);
end
fprintf('=====反幂法不使用Rayleigh quotient=====\n')
fprintf('特征值为： %0.7f\n',1.0/mu1+s)
fprintf('特征向量为 v = %0.7f,%0.7f,%0.7f\n',u0(1),u0(2),u0(3))
fprintf('迭代次数为： %d次\n',total)
clear all
A=[6,2,1;2,3,1;1,1,1];
s=6;
n=length(A(:,1));
u0=ones(n,1);%初始向量
v0=zeros(n,1);
%=====LU=求解Uv0=u0===
[B,p]=LUp(A-s*eye(n));
v0(n)=u0(n)/B(n,n);
for i=n-1:-1:1
    tmp=0;
    for j=i+1:n
        tmp=tmp+B(i,j)*v0(j);
    end
    v0(i)=(u0(i)-tmp)/B(i,i);
end
tol=1;
mu1=(v0'*u0)/sqrt(sum(u0.^2));
u0=v0/mu1;
total =1;
while(tol >0.0000001)
    a=mu1;
    tmp3=u0;
    %====Pu0 换行=======
    for i=1:n
        if(p(i)==i)
            continue;
        end
            tmp1=u0(i);
            u0(i)=u0(p(i));
            u0(p(i))=tmp1;
    end
    %===Ly=Pu0============
    for i=2:n
        tmp=0;
        for k=1:i-1
            tmp=tmp+B(i,k)*u0(k);
        end
        u0(i)=u0(i)-tmp;
    end
    %====Ux=y===========
    v0(n)=u0(n)/B(n,n);
    for i=n-1:-1:1
        tmp=0;
        for j=i+1:n
            tmp=tmp+B(i,j)*v0(j);
        end
        v0(i)=(u0(i)-tmp)/B(i,i);
    end
   mu1=(tmp3'*v0)/sqrt(sum(tmp3.^2));
   u0=v0/mu1;
   total =total +1;
   tol=abs(mu1-a);
end
fprintf('=====反幂法使用Rayleigh quotient=====\n')
fprintf('特征值为： %0.7f\n',1.0/mu1+s)
fprintf('特征向量为 v = %0.7f,%0.7f,%0.7f\n',u0(1),u0(2),u0(3))
fprintf('迭代次数为： %d次\n',total)
