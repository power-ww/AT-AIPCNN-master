function [SAM_SCC, SAM_groups,CC_groups,SAM_mat,SCC_index ] = SAM_SCC_zj( I_MS,I_HS,ratio,GNyq,flag_histmatch )
%   SAM_CC band assignment algorithm
%   Grouping of HS bands to fuse with MS imagery;
%   Input images are:
%   I_MS: MS images
%   I_HS: HS images
%   ratio: scale ratio between MS and HS images
%   GNyq: Gain at Nyquist frequency  GNyq = 0.29 .* ones(1,size(I_MS,3));
%   flag_histmatch: if 1, histogram matching is performed (default=1)

[L1,L2,Nb]=size(I_HS);
Nb_MS=size(I_MS,3);
cd functions/BandAssignment/Quality_indices_HS
I_PAN_LR=MTF_filter(I_MS,GNyq,ratio,41);
cd  ../Assignment
I_PAN_LR=I_PAN_LR(floor(ratio/2)+1:ratio:end,floor(ratio/2)+1:ratio:end,:);


if flag_histmatch~=0
    mean2_I_MS_LR=squeeze(sum(sum(I_HS,1),2))/L1/L2;
    std2_I_MS_LR=sqrt(squeeze(sum(sum(I_HS.^2,1),2))/L1/L2-mean2_I_MS_LR.^2);
    mean2_I_PAN_LR=squeeze(sum(sum(I_PAN_LR,1),2))/L1/L2;
    std2_I_PAN_LR=sqrt(squeeze(sum(sum(I_PAN_LR.^2,1),2))/L1/L2-mean2_I_PAN_LR.^2);
    I_MS_HM=zeros(L1,L2,Nb,Nb_MS);
    for ii1=1:Nb_MS
        for ii2=1:Nb
            I_MS_HM(:,:,ii2,ii1)=(I_HS(:,:,ii2)-mean2_I_MS_LR(ii2)*ones(L1,L2))*std2_I_PAN_LR(ii1)/std2_I_MS_LR(ii2)+mean2_I_PAN_LR(ii1)*ones(L1,L2);
        end
    end
else
    I_MS_HM=I_HS;
end

cd  ../Quality_Indices_HS
SCC_index=SCC_HSMS(  I_PAN_LR, I_MS_HM, flag_histmatch );
cd  ../Assignment
[~,SCC_max]=max(SCC_index);
CC_groups=cell(1,Nb_MS);
[M1,M2,Nb_MS]=size(I_MS);
[L1,L2,Nb]=size(I_HS);

if Nb_MS==1
    SAM_groups={1:Nb};
    return;
end

cd  ../Quality_Indices_HS
I_PAN_LR=MTF_filter(I_MS,GNyq,ratio,41);
I_PAN_LR=I_PAN_LR(floor(ratio/2)+1:ratio:end,floor(ratio/2)+1:ratio:end,:);
mean2_I_PAN_LR=squeeze(sum(sum(I_PAN_LR,1),2))/L1/L2;
std2_I_PAN_LR=sqrt(squeeze(sum(sum(I_PAN_LR.^2,1),2))/L1/L2-mean2_I_PAN_LR.^2);
cd  ../Assignment

mean2_I_MS_LR=squeeze(sum(sum(I_HS,1),2))/L1/L2;
std2_I_MS_LR=sqrt(squeeze(sum(sum(I_HS.^2,1),2))/L1/L2-mean2_I_MS_LR.^2);
I_PAN_HM=zeros(L1,L2,Nb_MS,Nb);
for ii1=1:Nb_MS
    for ii2=1:Nb
        I_PAN_HM(:,:,ii1,ii2)=(I_PAN_LR(:,:,ii1)-mean2_I_PAN_LR(ii1)*ones(L1,L2))*std2_I_MS_LR(ii2)/std2_I_PAN_LR(ii1)+mean2_I_MS_LR(ii2)*ones(L1,L2);
    end
end

cd  ../Quality_indices_HS
SAM_mat=SAM_HSMS(I_HS,I_PAN_HM);
[~,min_index]=min(SAM_mat);

a=ones(size(SCC_index));
b=a-SCC_index;
c=SAM_mat+b;
[~,SAMSCC_index]=min(c); 

cd  ../Assignment
SAM_groups=cell(1,Nb_MS);
SAM_SCC=cell(1,Nb_MS);
for ii1=1:Nb_MS
    CC_groups{ii1}=find(SCC_max==ii1);
    SAM_groups{ii1}=find(min_index==ii1);
    SAM_SCC{ii1}=find(SAMSCC_index==ii1);
end
cd ../../../
end

