function s=segcorr2_0319001(a,b,Y)
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

X1 = s1;
minv1 = min(X1(:)); maxv1 = max(X1(:));
X1 = (X1 - minv1) / (maxv1 - minv1);
X2 = s2;
minv2 = min(X2(:)); maxv2 = max(X2(:));
X2 = (X2 - minv2) / (maxv2 - minv2);
s=corr(X1(:),X2(:));