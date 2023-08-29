function X = recFUNC(dec,sINI,Lo,Hi,perFLAG)

% Reconstruction.
perm = [2,1,3];
W = cell(1,2);
for i = 1:2
    W{i} = wrec1D(dec{i,1},Lo{2},perm,perFLAG) + ...
        wrec1D(dec{i,2},Hi{2},perm,perFLAG);
end
X = (wrec1D(W{1},Lo{1},[],perFLAG) + wrec1D(W{2},Hi{1},[],perFLAG))/4;

% Extraction of central part
sREC = size(X);
F = floor((sREC-sINI)/2);
C = ceil((sREC-sINI)/2);
X = X(1+F(1):end-C(1),1+F(2):end-C(2),:);