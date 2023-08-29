% AT-AIPCNN for Hyperspectral image and multispectral image fusion
% Copyright(c) 2023 Xinyu Xu，Xiaojun Li
% All Rights Reserved.
%
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is here
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
%
% This is an implementation of the algorithm for A novel adatively
% optimized PCNN model for hyperspectral image sharpening.
%
% if you use this code, Please cite the following paper:
%
%  Xu, X.; Li, X.; Li, Y.; Kang, L.; Ge, J. A Novel Adaptively Optimized
%  PCNN Model for Hyperspectral Image Sharpening. Remote Sens. 2023, 15,
%  4205. https://doi.org/10.3390/rs15174205

clear
clc

load('data/d.mat');
sensor='none';
I_HS= imresize(I_HSS,ratio);
GNyq = 0.29 .* ones(1,size(I_MS,3));
L=11;Q_blocks_size=32;flag_cut_bounds=1;dim_cut=11;th_values=0;bicubic=0;
[M,N,H]=size(I_MS);

addpath(genpath('AT-AIPCNN'))
addpath(genpath('functions'))
addpath(genpath('functions/CSA'))
addpath(genpath('functions/ONR'))
addpath(genpath('functions/BandAssignment'))
addpath(genpath('functions/Quality_Indices'))

I_MS_norm=norm21(I_MS);
I_HS_LR_norm=norm21(I_HSS);
[SAM_SCC, SAM_groups,CC_groups,SAM_mat,SCC_index ] = SAM_SCC_zj(I_MS_norm,I_HS_LR_norm,ratio,GNyq,1 );

Z=length(SAM_SCC);
MS=cell(1,Z);
HS=cell(1,Z);
I_HS_LR=cell(1,Z);
I_GT=cell(1,Z);
for i=1:Z
    MS{i}=I_MS(:,:,i);
    HS{i}=I_HS(:,:,SAM_SCC{1,i});
    I_HS_LR{i}=I_HSS(:,:,SAM_SCC{1,i});
    I_GT{i}=REF(:,:,SAM_SCC{1,i});
end

Z=length(HS);
HSHS=cell(1,Z);
GT=cell(1,Z);
for i=1:Z
    if length(SAM_SCC{i})>5
        A{i} = HS{i};
        B{i} = I_GT{i};
        
        X = permute(A{i}, [3, 1, 2]);  
        X = X(:, :);
        minv = min(X(:)); maxv = max(X(:));
        X = (X - minv) / (maxv - minv);

        k = 5; 
        ONR_L = ONR_init(X');
        band_set = ONR(X, ONR_L, k);
        HSHS{i}=A{i}(:,:,band_set);
        GT{i}=B{i}(:,:,band_set);
    end
    if isempty(HSHS{i})
        HSHS{i}=HS{i};
        GT{i}=I_GT{i};
    end
end

Function_name='F10';
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);
pop_num=20;
Max_iter=30;
fMin_CSA=zeros(Z,1);
bestX_CSA=zeros(Z,dim);
CSA_curve=zeros(Z,Max_iter);
bestX_SSA=zeros(Z,dim);

t1=clock;
for i=1:Z
    [fMin_CSA(i),bestX_CSA(i,:),CSA_curve(i,:)]=CSA_zj(pop_num,Max_iter,lb,ub,dim,fobj,HSHS{i},MS{i},ratio,GT{i}); %%F7 PCNNATWT_0307
end
t2=clock;
time_CSA=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);
for i=1:Z
    bestX_SSA(i,:)=bestX_CSA(i,:);
end

thresholdvalue=0;
I_NSCT_AWLP=cell(1,Z);
t1=clock;
for i=1:Z
    [~,I_NSCT_AWLP{i} ] = PCNNATWT_0329_002(HS{i},MS{i},ratio,bestX_SSA(i,1),bestX_SSA(i,2),bestX_SSA(i,3),bestX_SSA(i,4),bestX_SSA(i,5));
end
t2=clock;
time_fusion=(t2(end)+t2(end-1)*60+t2(end-2)*3600-t1(end)-t1(end-1)*60-t1(end-2)*3600);

I_ATAIPCNN=zeros(size(I_HS));
for i=1:Z
    I_ATAIPCNN(:,:,SAM_SCC{1,i})=I_NSCT_AWLP{1,i};
end
I_ATAIPCNN=I_ATAIPCNN./max(I_ATAIPCNN(:));
REF=REF./max(REF(:));

[psnr, rmse, ergas, sam, uiqi, ssim, DD, CCS] = quality_assessment(double(im2uint8(I_ATAIPCNN)),double(im2uint8(REF)),0, 1/ratio);
Results(1,:)= [psnr,rmse, ergas, sam, uiqi,ssim,DD,CCS] ;