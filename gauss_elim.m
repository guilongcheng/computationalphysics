function[A,b,x]=gauss_elim(A,b)
n=length(b);
x=zeros(1,n);

for i=1:n-1
    if(A(i,i)==0)
        t=min(find(A(i+1:n)~=0)+i);
        if(isempty(t))
            disp("gauss_elim error:A matrix is singular");
            return
        end;
        temp=A(i,:);tb=b(i);
        A(i,:)=A(t,:);b(i)=b(t);
        A(t,:)=temp;b(t)=tb;
    end;
    for j=i+1:n
        m=-A(j,i)/A(i,i);
        A(j,i)=0;
        A(j,i+1:n)=A(j,i+1:n) + m*A(i,i+1:n);
        b(j)=b(j)+m*b(i);
    end;
end;

x(n)=b(n)/A(n,n);
for i = n-1:-1:1
    x(i)=(b(i)-sum(x(i+1:n).*A(i,i+1:n)))/A(i,i);
end;

