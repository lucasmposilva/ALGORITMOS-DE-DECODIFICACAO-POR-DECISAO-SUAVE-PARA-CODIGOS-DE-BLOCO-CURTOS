function INPUT_SD_CHASE_II_BM_BCH(N, K, dmin, t, radius, input_file, output_method)
% generator polynomial
[gpol, t] = bchgenpoly(N,K,[],'double');

% generating the set of test patterns
tps = chaseIITpsGen(radius);

% BCH enc and decod
enc = comm.BCHEncoder(N,K,gpol);
dec = comm.BCHDecoder(N,K,gpol);

% matrix to count words decoded by HD
nSkipped = zeros(1,1);

% BPSK mod and demod
modulator = comm.BPSKModulator;
demod = comm.BPSKDemodulator;

% results matrices
load('results_header.mat');
results = -1*ones(1,size(header,2));

% loading data var from input file
load(input_file,'data');
data = data';

% adjusting parameters to parse input data
nBits    = size(data,1).*size(data,2);
tailBits = mod(size(data,1).*size(data,2),N);
nBlocks  = ceil(size(data,1).*size(data,2) ./ N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error rates
channelErrorRate = comm.ErrorRate('ComputationDelay',0);
decodingErrorRate = comm.ErrorRate('ComputationDelay',0);
% Counters
n=0;ne=0;

tic
% for each loop(msg block)
for blockMsgKey = 1:nBlocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Source
    msg        = zeros(K, 1);
    encodedMsg = zeros(N, 1);
    
    % Rx
    if(blockMsgKey==nBlocks)
        r   = [ data( ((blockMsgKey-1)*N)+1 : ((blockMsgKey-1)*N)+tailBits)'; ...
                ones(N-tailBits, 1)];
    else    
        r   = [ data( ((blockMsgKey-1)*N)+1 : blockMsgKey*N )' ];
    end
    rHD = step(demod, r);
    r   = real(r);
        
    % Decoding
    [guess, x] = step(dec, rHD);
    ccHD = step(enc, guess);

    sumAlt = 0; sumRef = 0;

    if ccHD ~= -1
        altBits = xor(rHD,ccHD);
    else
        altBits = zeros(N,1);
    end
    nAltBits = sum(altBits);
    sumAlt = altBits' * abs(r);

    if(dmin - nAltBits < radius)
        diff = dmin - nAltBits;
    else
        diff = radius;
    end
    values = zeros(1, radius);
    indexes = zeros(1, radius);
    % grant that the altered bits won't be changed by the tps
    rTemp = (10*altBits)+r;
    % looking for the "radius" lowest values
    for x = 1:radius
        [values(x), indexes(x)] = min(abs(rTemp));
        rTemp(indexes(x)) = 10;
    end
    sumRef = sum(values(1:diff));
    commitedCand = 0;

    if ((sumRef>=sumAlt) && (size(ccHD,1)==N))
        guess = ccHD;
        nSkipped(1) = nSkipped(1) + 1;
    else            
        % building candidates
        nCand = size(tps,2);
        cand = repmat(rHD,1,nCand);
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

        % Correlation
        cand = 1 - 2*cand(:,1:commitedCand);
        Correlation = r' * cand(:,1:commitedCand);    
        [MaxC, GreaterC] = max(Correlation);
        guess = cand(:,GreaterC);
        guess = (-guess+1)./2;
    end

    % Errors computation
    channelErrorStats = step(channelErrorRate, encodedMsg, rHD);
    decodingErrorStats = step(decodingErrorRate, encodedMsg, guess);
    if( sum( xor(encodedMsg,guess) ) )
        ne = ne+1;
    end
    n = n + 1;        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOP end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% logging loop time
deltaT = toc;
disp("Simulation time: "+deltaT+" s")

% writing into results variable
results(1,1:12) = [ ...
    -1, deltaT, ...
    n, nSkipped(1), n-nSkipped(1), ne, ne/n, ...
    decodingErrorStats(3), ...
    decodingErrorStats(2), decodingErrorStats(1), ...
    channelErrorStats(2), channelErrorStats(1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saving results
saveResults(output_method,results,header,["SD_CHASE_II_BM_BCH", input_file]);