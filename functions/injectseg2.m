function FusedBlock=injectseg2(I_MS,DetailsHRPan,b,LG,Y)

[m,n]=size(I_MS(:,:,b));
FusedBlock=zeros(m,n);

for i=1:m
    for j=1:n
         if Y(i,j)==1 
             
            BlockDetails = DetailsHRPan(i,j);
            Details2Add=LG*BlockDetails;
            FusedBlock(i,j) = I_MS(i,j,b) + Details2Add;      
        end
    end
end

