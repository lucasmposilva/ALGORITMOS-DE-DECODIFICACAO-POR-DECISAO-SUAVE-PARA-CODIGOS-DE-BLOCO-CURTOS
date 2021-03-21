function [cMatrix, commited] = insertCandNoCheck(cMatrix,cand,commited)
%
cMatrix(:,commited+1) = cand;
commited = commited + 1;
