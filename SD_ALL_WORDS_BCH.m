function SD_ALL_WORDS_BCH_15_7(vEbNo, wethr, output_method)
% parameters of the code
N = 15; K = 7; dmin = 5; t = 2;

% loading g, h and s matrices
load('env_BCH_15_7.mat');

% getting number of rounds, one index per round 
Nindexes = size(vEbNo, 2);

% BPSK mod and demod
mod   = comm.BPSKModulator;
demod = comm.BPSKDemodulator;

% loading all words
load('words_bch_15_7.mat');

% BPSK modulation
words = 1 - 2*words;

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
        encodedData = encode(data,N,K,'linear/binary',g);
        modSignal = step(mod, encodedData);
        r = step(chan, modSignal);
        
        % Rx
        rHD = step(demod, r);
        r = real(r);
        
        % Decoding
        Correlation = r' * words;    
        [MaxC, GreaterC] = max(Correlation);
        
        guess = words(:,GreaterC);
        guess = (-guess+1)./2;

        decodingErrorStats = step(decodingErrorRate, encodedData, guess);
        channelErrorStats = step(channelErrorRate, encodedData, rHD);
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
    disp("EbNo: "+EbNo+"dB"+" Simulation time: "+deltaT+" s");
    
    % writing into results variable
    results(INDEX,[1:12]) = [ ...
        vEbNo(INDEX), deltaT, n, 0, n,...
        ne, ne/n, ...
        decodingErrorStats(3), ...
        decodingErrorStats(2), decodingErrorStats(1) ...
        channelErrorStats(2), channelErrorStats(1)];
% ROUND end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saving results
saveResults(output_method,results,header,["SD_ALL_WORDS_BCH_15_7",mat2str(vEbNo)]);