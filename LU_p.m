function [A,p]=LUp(A)
    %选主元的三角分解方法，B储存LU，p为换行信息
    %matlab函数传值而不传地址，因此不会改变A
    n=length(A(:,1));
    p=zeros(n,1);
    for k=1:n
       for i=k:n
    	   A(i,k)=A(i,k)-A(i,1:k-1)*A(1:k-1,k);
       end
       [m,p(k)]=max(abs(A(k:n,k)));
       p(k)=p(k)+k-1;
       if(p(k)~=k)
           tmp=A(k,:);
           A(k,:)=A(p(k),:);
           A(p(k),:)=tmp;
       end
       if(k~=n)
            A(k+1:n,k)=A(k+1:n,k)/A(k,k);
       end
       if(k~=n)
           for i=k+1:n
                 A(k,i)=A(k,i)-A(k,1:k-1)*A(1:k-1,i);
           end
       end
    end
end
