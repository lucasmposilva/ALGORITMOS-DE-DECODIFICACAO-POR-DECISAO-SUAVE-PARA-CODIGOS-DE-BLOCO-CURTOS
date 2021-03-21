function INPUT_SD_ALL_WORDS_HAMMING(input_file, output_method)
% parameters of the code
N = 7; K = 4; dmin = 3; t = 1;

% loading g, h and s matrices
load('env_HAMMING.mat')

% BPSK mod and demod
% mod   = comm.BPSKModulator;
demod = comm.BPSKDemodulator;

% loading all words
load('words_hamming.mat')

% BPSK modulation
words = 1 - 2*words;

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
    Correlation = r' * words;    
    [MaxC, GreaterC] = max(Correlation);

    guess = words(:,GreaterC);
    guess = (-guess+1)./2;

    decodingErrorStats = step(decodingErrorRate, encodedMsg, guess);
    channelErrorStats = step(channelErrorRate, encodedMsg, rHD);
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
disp("Simulation time: "+deltaT+" s");

% writing into results variable
results(1,1:12) = [ ...
    -1, deltaT, ...
    n, 0, n, ne, ne/n, ...
    decodingErrorStats(3), ...
    decodingErrorStats(2), decodingErrorStats(1), ...
    channelErrorStats(2), channelErrorStats(1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saving results
saveResults(output_method,results,header,["SD_ALL_WORDS_HAMMING_INPUT",input_file]);