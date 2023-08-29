% AT-AIPCNN for Hyperspectral image and multispectral image fusion
% Copyright(c) 2023 Xinyu Xu,Xiaojun Li
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

function [Correlationci,I_Fus_AWLP] = PCNNATWT_0329_002(I_MS,I_PAN,ratio,al,ae,af,B,w)

[Height,Width,Bands]=size(I_MS);
I_Fus_AWLP=zeros(Height,Width,Bands,'double');

Correlationci=[];

I_PAN = repmat(I_PAN,[1 1 size(I_MS,3)]);

for ii = 1 : size(I_MS,3)
    I_PAN(:,:,ii) = (I_PAN(:,:,ii) - mean2(I_PAN(:,:,ii))).*(std2(I_MS(:,:,ii))./std2(I_PAN(:,:,ii))) + mean2(I_MS(:,:,ii));
end

h=[1 4 6 4 1 ]/16;
g=[0 0 1 0 0 ]-h;
htilde=[ 1 4 6 4 1]/16;
gtilde=[ 0 0 1 0 0 ]+htilde;
h=sqrt(2)*h;
g=sqrt(2)*g;
htilde=sqrt(2)*htilde;
gtilde=sqrt(2)*gtilde;
WF={h,g,htilde,gtilde};
Levels = ceil(log2(ratio));

for b=1:Bands
    
    WT = ndwt2_working(I_PAN(:,:,b),Levels,WF);
    for ii = 2 : numel(WT.dec)
        WT.dec{ii} = zeros(size(WT.dec{ii}));
    end
    
    imageHRLP =newindwt2(WT,'c');
    StepDetails = I_PAN(:,:,b) - imageHRLP;
    
    [m,n]=size(I_MS(:,:,b));
    
    F=EstimateSpatialFrequency(I_MS(:,:,b),5);
    F=norm21(F);
    
    I1=EstimateSignalVarianceLocally(I_MS(:,:,b),5);
    I1=(norm21(I1)+norm21(I_MS(:,:,b)))/2;
    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    sum1=0;
    mt=1;
    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    %^^^^^^^^^^^^^^^^^^     Initialize    PCNN     Parameters     ^^^^^^^^^^^^^
    % P:    ¦ÁL	     ¦ÁE	     ¦ÁF	       £ÖF	    £ÖL  	   £ÖE	  ¦Â
    % V:    1£®0	     1£®0     0£®1	   0£®5	    0£®2	       20	  0£®1
    %        al=1.0;  ae=0.62;   af=0.1;    vf=0.5;   vl=0.2;     ve=20; B=0.1;
    %        ae=0.62;
    vf=0.5;    vl=0.2;     ve=20;
    %^^^^^^^^^^^^^^^^^^
    W=[w 1 w;
        1 0 1;
        w 1 w];
    % W=[0.707 1 0.707;
    %     1  0  1;
    %    0.707 1 0.707];
    %^^^^^^^^^^^^^^^^^^
    M=W;     Y=zeros(m,n);
    L=Y;       U=Y;
    E=ve*ones(m,n);
    %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    while mt
        
        %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        %^^^^^^^^^^^^^^^^^^           test     PCNN           ^^^^^^^^^^^^^
        K=conv2(Y,M,'same');
        F=exp(-af)*F+vf*K+I1;
        L=exp(-al)*L+vl*K;
        U=F.*(1+B*L);
        Y=double(U>E);
        E=exp(-ae)*E+65535*Y;
        
        sum1=sum1+sum(Y(:));
        if sum1>=(m*n)
            mt=0;
        end
        %^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        YY=ones(Height,Width);
        StdMS_all=segstd2(I_MS(:,:,b),YY);
        StdPan_all=segstd2(imageHRLP,YY);
        all=StdMS_all/StdPan_all;
        if sum(Y(:))~=0
            StdMS=segstd2(I_MS(:,:,b),Y);
            StdPan=segstd2(imageHRLP,Y);
            
            Corr=segcorr2_0319001(I_MS(:,:,b),imageHRLP,Y);
            
            Correlationtemp1=segcov2(I_MS(:,:,b),I_PAN(:,:,b),Y);
            Correlationtemp2=segcov2(imageHRLP,I_PAN(:,:,b),Y);
            Correlation=Correlationtemp1(1,2)/(Correlationtemp2(1,2)+eps);
            if (Correlation>0)
                LG=Corr*((StdMS/StdPan)/all);
            else
                LG=1;
            end
            
            temp1=injectseg2(I_MS,StepDetails,b,LG,Y);
            temp2=I_Fus_AWLP(:,:,b)+temp1;
            I_Fus_AWLP(:,:,b)=temp2;
        end
        
    end
    
end

end