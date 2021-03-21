function [tps] = chaseIITpsGen(r)
%

tps = [de2bi([0:(2.^r)-1])]';