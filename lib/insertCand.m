function [cMatrix, commited] = insertCand(cMatrix,cand,commited)
%
ans = (size(cand,1)==sum(cand==cMatrix));
if ( (sum(ans(1,:))==0) || commited==0 )
    cMatrix(:,commited+1) = cand;
    commited = commited + 1;
end