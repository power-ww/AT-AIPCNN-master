function X = wrec1D(X,F,perm,perFLAG)

if ~isempty(perm) , X = permute(X,perm); end
if perFLAG
    nb = length(F)-1;
    X = [X X(:,1:nb,:)];
end
X = convn(X,F);
if ~isempty(perm) , X = permute(X,perm); end