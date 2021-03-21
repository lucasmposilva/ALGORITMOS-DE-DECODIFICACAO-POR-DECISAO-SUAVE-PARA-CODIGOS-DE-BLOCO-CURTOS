function bchWordsGen(n,k,name)
    
    [gpol, t]   = bchgenpoly(n,k,[],'double');
    enc         = comm.BCHEncoder(n,k,gpol);
    
    words = zeros(n, 2^(k));
    
    vMsg  = de2bi([0:2^(k)-1])';
    
    for i = 1:size(vMsg,2)
        words(:,i)  = step(enc, vMsg(:,i));
    end
    
    save(name,'words');