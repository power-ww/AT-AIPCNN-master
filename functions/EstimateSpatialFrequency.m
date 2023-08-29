function FFF= EstimateSpatialFrequency(Y,m)
n=m;
x1=Y;
[height, width]=size(Y);
FFF=zeros(height, width);
for i=1:height-n+1
    for j=1:width-n+1
        c=x1(i:i+(n-1),j:j+(n-1));
        I2=c;M=n;N=n;
        
        SumRF=0;
        for ii=1:M
            for jj=2:N
                SumRF = SumRF + (I2(ii,jj)-I2(ii,jj-1))^2;
            end
        end
        RF=sqrt(SumRF/(M*N));

        SumCF=0;
        for ii=2:M
            for jj=1:N
                SumCF = SumCF + (I2(ii,jj)-I2(ii-1,jj))^2;
            end
        end
        CF=sqrt(SumCF/(M*N));
        
        SF=sqrt(RF^2+CF^2);
        FFF(i+(n-1)/2,j+(n-1)/2)=SF/(n*n);
    end
end
end