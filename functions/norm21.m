function s=norm21(a)

amax=max(a(:));
amin=min(a(:));

s=(a-amin)./(amax-amin);
