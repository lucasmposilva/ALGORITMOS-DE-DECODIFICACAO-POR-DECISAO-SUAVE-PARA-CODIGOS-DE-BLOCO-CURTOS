function SD_CHASE_II_BM_BCH(N, K, dmin, t, radius, vEbNo, wethr, output_method)
% generator polynomial
[gpol, t] = bchgenpoly(N,K,[],'double');

% generating the set of test patterns
tps = chaseIITpsGen(radius);

% BCH enc and decod
enc = comm.BCHEncoder(N,K,gpol);
dec = comm.BCHDecoder(N,K,gpol);

% getting number of rounds, one index per round 
Nindexes = size(vEbNo, 2);

% matrix to count words decoded by HD
nSkipped = zeros(1,Nindexes);

% BPSK mod and demod
mod = comm.BPSKModulator;
demod = comm.BPSKDemodulator;

% results matrices
load('results_header.mat');
results = -1*ones(Nindexes, size(header,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each round(EbNo)
for INDEX = 1:1:Nindexes
    EbNo = vEbNo(INDEX);
    % Adjusting for correct power distribution
    EbNoEff = EbNo + 10*log10(K/N);
    % Channel parameters
    chan = comm.AWGNChannel(...
        'NoiseMethod','Signal to noise ratio (Eb/No)','EbNo',EbNoEff);
    % Error rates
    channelErrorRate = comm.ErrorRate('ComputationDelay',0);
    decodingErrorRate = comm.ErrorRate('ComputationDelay',0);
    % Counters
    n=0;ne=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    % for each loop(msg)
    while ne < wethr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Source
        data = randi([0 1], K, 1);
        
        % Tx
        encodedData = step(enc, data);
        modSignal = step(mod, encodedData);
        r = step(chan, modSignal);
        
        % Rx
        rHD = step(demod, r);
        r = real(r);
        
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
            nSkipped(INDEX) = nSkipped(INDEX) + 1;
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
        channelErrorStats = step(channelErrorRate, encodedData, rHD);
        decodingErrorStats = step(decodingErrorRate, encodedData, guess);
        if( sum( xor(encodedData,guess) ) )
            ne = ne+1;
        end
        n = n + 1;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOOP end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % logging loop time
    deltaT = toc;
    disp("EbNo: "+EbNo+"dB"+" Simulation time: "+deltaT+" s")
    
    % writing into results variable
    results(INDEX,1:12) = [ ...
        vEbNo(INDEX), deltaT, ...
        n, nSkipped(INDEX), n-nSkipped(INDEX), ...
        ne, ne/n, ...
        decodingErrorStats(3), ...
        decodingErrorStats(2), decodingErrorStats(1), ...
        channelErrorStats(2), channelErrorStats(1), ...
        ];
% ROUND end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saving results
saveResults(output_method,results,header,["SD_CHASE_II_BM_BCH",mat2str(vEbNo)]);