function [lb,ub,dim,fobj] = Get_Functions_details(F)
switch F
        case 'F10'       
        fobj = @F10;
        lb=[0.1,0.1,0.1,0.1,0.1];
        ub=[3,3,3,1,1];
        dim=5; 
end
end

% F10
function o = F10(X,I_MS,I_PAN,ratio,I_GT)
al=X(1);
ae=X(2);
af=X(3);
B=X(4);
w=X(5);

[~,I_F] = PCNNATWT_0329_002(I_MS,I_PAN,ratio,al,ae,af,B,w);
I_GT=I_GT./max(I_GT(:));
I_F=I_F./max(I_F(:));
x=I_GT;
y=I_F;
if size(I_GT,3)>1
    sz_x = size(x);
    n_bands = sz_x(3);
    n_samples = sz_x(1)*sz_x(2);
else
    sz_x = size(x);
    n_bands = 1;
    n_samples = sz_x(1)*sz_x(2);
end
aux = sum(sum((x - y).^2), 2)/n_samples;
rmse_per_band = sqrt(aux);

sam= SpectAngMapper(double(im2uint8(I_F)),double(im2uint8(I_GT)));
sam_max=round(max(sam));
sam_min=round(min(sam));
sam_change=sam_max-sam_min;

mean_y = sum(sum(y, 1), 2)/n_samples;
ratio_ergas=1/9;
ergas = 100*ratio_ergas*sqrt(sum((rmse_per_band ./ mean_y).^2)/n_bands);
ergas_max=round(max(ergas));
ergas_min=round(min(ergas));
ergas_change=ergas_max-ergas_min;

to=sam_change+ergas_change;
a=sam_change/(to+eps);
b=ergas_change/(to+eps);
if a==0 || b==0
    o=sam+ergas;
else
    o=b*sam+a*ergas;
end

end
