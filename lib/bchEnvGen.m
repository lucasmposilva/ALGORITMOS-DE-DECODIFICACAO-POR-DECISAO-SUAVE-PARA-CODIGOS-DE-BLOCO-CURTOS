function [g, h, s, t] = bchEnvGen(n, k, pp)

if(nargin > 3 || nargin < 2 )
    disp("See the funtion help")
else
    if(nargin == 2)
        [gpol, t] = bchgenpoly(n,k);
    elseif(nargin == 3)
        [gpol, t] = bchgenpoly(n,k,pp);
    end  
    [g, h]  = cyclgen(n,double(gpol.x));
    s       = syndtable(h);
end