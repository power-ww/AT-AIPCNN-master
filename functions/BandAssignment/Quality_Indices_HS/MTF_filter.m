function [ I_out ] = MTF_filter( I_in,MTF_Nyq,scale,N )

if nargin<4 
    N=scale*10+1;
end
I_in=double(I_in);

[L1,L2,Nb]=size(I_in);
I_out=zeros(L1,L2,Nb);
for ii=1:Nb
    sigma=sqrt(-((N-1)/scale/2)^2/2/log(MTF_Nyq(ii)));
    filter=fspecial('gaussian',N,sigma);
    filter=fwind1(filter./max(filter(:)),kaiser(N));
    I_out(:,:,ii)=imfilter(I_in(:,:,ii),real(filter),'replicate');
end
end

