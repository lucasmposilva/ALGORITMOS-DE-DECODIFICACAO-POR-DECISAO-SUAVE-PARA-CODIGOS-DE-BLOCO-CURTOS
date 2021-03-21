function [w, cc] = myDec(cw,h,s,n,k)
%
cw = cw';
syn = cw*h';
syn = mod(syn,2);
Nsyn = bi2de(flip(syn)) + 1;
%
if Nsyn<2^(n-k) || Nsyn >=0
    cc = mod(cw + s(Nsyn,:),2)';
    w = cc(n-k+1:n);
else
    cc = -1;
    w = cc;
end
end
%