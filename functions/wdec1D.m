function [L,H] = wdec1D(X,Lo,Hi,perm,dwtEXTM)

if ~isempty(perm) , X = permute(X,perm); end
sX = size(X);
if length(sX)<3 , sX(3) = 1; end
lf = length(Lo);
lx = sX(2);
lc = lx+lf-1;
switch dwtEXTM
    case 'zpd'             % Zero extension.
        
    case {'sym','symh'}    % Symmetric extension (half-point).
        X = [X(:,lf-1:-1:1,:) , X , X(:,end:-1:end-lf+1,:)];
        
    case 'sp0'             % Smooth extension of order 0.
        X = [X(:,ones(1,lf-1),:) , X , X(:,lx*ones(1,lf-1),:)];
        
    case {'sp1','spd'}     % Smooth extension of order 1.
        Z = zeros(sX(1),sX(2)+ 2*lf-2,sX(3));
        Z(:,lf:lf+lx-1,:) = X;
        last = sX(2)+lf-1;
        for k = 1:lf-1
            Z(:,last+k,:) = 2*Z(:,last+k-1,:)- Z(:,last+k-2,:);
            Z(:,lf-k,:)   = 2*Z(:,lf-k+1,:)- Z(:,lf-k+2,:);
        end
        X = Z; clear Z;
        
    case 'symw'            % Symmetric extension (whole-point).
        X = [X(:,lf:-1:2,:) , X , X(:,end-1:-1:end-lf,:)];
        
    case {'asym','asymh'}  % Antisymmetric extension (half-point).
        X = [-X(:,lf-1:-1:1,:) , X , -X(:,end:-1:end-lf+1,:)];        
        
    case 'asymw'           % Antisymmetric extension (whole-point).
        X = [-X(:,lf:-1:2,:) , X , -X(:,end-1:-1:end-lf,:)];

    case 'rndu'            % Uniformly randomized extension.
        X = [randn(sX(1),lf-1,sX(3)) , X , randn(sX(1),lf-1,sX(3))];        
                        
    case 'rndn'            % Normally randomized extension.
        X = [randn(sX(1),lf-1,sX(3)) , X , randn(sX(1),lf-1,sX(3))];        
                
    case 'ppd'             % Periodized extension (1).
        X = [X(:,end-lf+2:end,:) , X , X(:,1:lf-1,:)];
        
    case 'per'             % Periodized extension (2).
        if rem(lx,2) , X = [X , X(:,end,:)]; end
        X = [X(:,end-lf+2:end,:) , X , X(:,1:lf-1,:)];        
end
L = convn(X,Lo);
H = convn(X,Hi);
clear X
switch dwtEXTM
    case 'zpd'
    otherwise
        lenL = size(L,2);
        first = lf; last = lenL-lf+1;
        L = L(:,first:last,:); H = H(:,first:last,:);
        lenL = size(L,2);
        first = 1+floor((lenL-lc)/2);  last = first+lc-1;
        L = L(:,first:last,:); H = H(:,first:last,:);
end
if isequal(dwtEXTM,'per')
    first = 1; last = lx;
    L = L(:,first:last,:);
    H = H(:,first:last,:);
end

if ~isempty(perm)
    L = permute(L,perm);
    H = permute(H,perm);
end