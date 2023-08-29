function [ CC_groups,SCC_index ] = CC_classifier_v2( I_MS,I_HS,ratio,GNyq,flag_histmatch )
%CC_CLASSIFIER Grouping HS images for fusing with MS
%   I_MS_LR: original HS images
%   I_PAN:   original MS images
%   ratio:   scale ratio between MS and HS images
%   GNyq:    MTF gain at Nyquist frequency
%   flag_histmatch: if 1, histogram matching is performed (default=0)

if nargin<=4
    flag_histmatch=0;
end

[L1,L2,Nb]=size(I_HS);
Nb_MS=size(I_MS,3);
cd ../Quality_indices_HS
I_PAN_LR=MTF_filter(I_MS,GNyq,ratio,41);
cd ../Assignment
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

cd ../Quality_Indices_HS
SCC_index=SCC_HSMS(  I_PAN_LR, I_MS_HM, flag_histmatch );
cd ../Assignment

[~,SCC_max]=max(SCC_index);

CC_groups=cell(1,Nb_MS);
for ii1=1:Nb_MS
    CC_groups{ii1}=find(SCC_max==ii1);
end

end

