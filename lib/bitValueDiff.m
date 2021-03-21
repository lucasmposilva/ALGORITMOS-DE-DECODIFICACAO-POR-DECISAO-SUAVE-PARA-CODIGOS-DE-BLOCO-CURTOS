function [ans] = bitValueDiff(mCand,Ncand)
%
ref  = sum(mCand,2);
vONES   = (ref==Ncand);
vZEROS  = (ref==0);
ans = sum(not(xor(vONES,vZEROS)));
%