function s=segstd2(a,Y)
[m,n]=size(a);
num=sum(Y(:));
k=1;
s1=zeros(1,num);
for i=1:m
    for j=1:n
        if Y(i,j)==1
            s1(1,k)=a(i,j);
            k=k+1;
        end
    end
end
X = s1;
minv = min(X(:)); maxv = max(X(:));
X = (X - minv) / (maxv - minv);
s=std(X);