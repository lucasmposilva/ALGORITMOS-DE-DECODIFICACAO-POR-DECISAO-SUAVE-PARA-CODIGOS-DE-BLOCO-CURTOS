function FER_BER_CH(K, vEbNo, wethr, output_method)
% getting number of rounds, one index per round 
Nindexes = size(vEbNo, 2);

% BPSK mod and demod
mod   = comm.BPSKModulator;
demod = comm.BPSKDemodulator;

% results matrices
load('results_header.mat');
results = -1*ones(Nindexes,size(header,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each round(EbNo)
for INDEX = 1:1:Nindexes
    EbNo = vEbNo(INDEX);
    % Adjusting for correct power distribution
    EbNoEff = EbNo;
    % Channel parameters
    chan = comm.AWGNChannel(...
        'NoiseMethod','Signal to noise ratio (Eb/No)','EbNo',EbNoEff);
    % Error rates
    channelErrorRate = comm.ErrorRate('ComputationDelay',0);
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
        encodedData = data;
        modSignal   = step(mod, encodedData);
        r           = step(chan, modSignal);
        
        % Rx
        rHD = step(demod, r);
        
        % Channel errors computation
        channelErrorStats = step(channelErrorRate, encodedData, rHD);
        %
        if( sum( xor(data,rHD) ) )
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
    results(INDEX,1:12) = [ ...
        vEbNo(INDEX), deltaT, ...
        n, n, 0, ne, ne/n, ...
        channelErrorStats(3), ...
        -1, -1, ...
        channelErrorStats(2), channelErrorStats(1)];
% ROUND end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saving results
saveResults(output_method,results,header,["BER_FER_CH",mat2str(vEbNo)]);