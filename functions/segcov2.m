function s=segcov2(a,b,Y)
[m,n]=size(a);
num=sum(Y(:));
k=1;
s1=zeros(1,num);
s2=zeros(1,num);
for i=1:m
    for j=1:n
        if Y(i,j)==1
            s1(1,k)=a(i,j);
            s2(1,k)=b(i,j);
            k=k+1;
        end
    end
end

s=cov(s1(:),s2(:));
if num==0
    s=zeros(2,2);
end