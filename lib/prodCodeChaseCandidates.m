function [mCand, map2b0, map2b1] = prodCodeChaseCandidates(nWords,n,k,gpol,radius,rMat)

% generating the set of test patterns
tps = chaseIITpsGen(radius);

% encoder
enc = comm.BCHEncoder(n,k,gpol);
% decoder
dec = comm.BCHDecoder(n,k,gpol);

for i = 1:nWords
    % The round word
    r = real( rMat(i,:)' );
    rHD = (r < 0);

    % Decoding
    [guess, x] = step(dec, rHD);
    ccHD = step(enc, guess);

    if ccHD ~= -1
        altBits = xor(rHD,ccHD);
    else
        altBits = zeros(N,1);
    end
    
    values = zeros(1, radius);
    indexes = zeros(1, radius);
    % grant that the altered bits won't be changed by the tps
    rTemp = r; %+ (10*altBits);
    % looking for the "radius" lowest values
    for x = 1:radius
        [values(x), indexes(x)] = min(abs(rTemp));
        rTemp(indexes(x)) = 10;
    end
    commitedCand = 0;
            
    % building candidates
    nCand = size(tps,2);
    cand = repmat(ccHD,1,nCand);
    % for each tp
    for x = 1:nCand
        temp = cand(:,x);
        temp(indexes) = tps(:,x);
        %
        [decDataTP, nc] = step(dec, temp);
        a = step(enc, decDataTP);
        % if a correction is possible
        if (a ~= -1)
            [cand, commitedCand] = insertCandNoCheck(cand,a,commitedCand);
        end
    % end for each tp
    end
    % Mapping
    mCand(:,:,i)    = 1 - 2*cand(:,:);
    map2b0(:,:,i)   = (cand == 0);
    map2b1(:,:,i)   = cand;
end