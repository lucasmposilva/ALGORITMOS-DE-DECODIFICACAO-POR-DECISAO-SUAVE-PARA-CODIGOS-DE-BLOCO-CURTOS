function INPUT_HD_BM_BCH(N,K,input_file, method, output_method)
%
% function INPUT_HD_BM_BCH(N, K, input_file, method, output_method)
%
% expected var from input_file: data
%
% N: word length
% K: message length
%
% method:
%   'column'      - parse data by column
%   'row'         - parse data by row
%
% output_method:
%   'a' - append 
%   'w' - write
%

% generator polynomial
[gpol, t] = bchgenpoly(N,K,[],'double');

% BPSK demod
demod = comm.BPSKDemodulator;

% BCH enc and decod
enc = comm.BCHEncoder(N,K,gpol);
dec = comm.BCHDecoder(N,K,gpol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error rates
channelErrorRate = comm.ErrorRate('ComputationDelay',0);
decodingErrorRate = comm.ErrorRate('ComputationDelay',0);
% Counters
n=0;ne=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
% loading data var from input file
load(input_file,'data');
% adjusting data variable
if(method=="row")
    data = data';
elseif(method~="column")
    disp('Invalid method...');
    quit;
end
% results matrices
load('results_header.mat');
results = -1*ones(1,size(header,2));
% adjusting parameters to parse input data
nBits    = size(data,1).*size(data,2);
tailBits = mod(size(data,1).*size(data,2),N);
nBlocks  = ceil(size(data,1).*size(data,2) ./ N);
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
    [guess, cnumerror] = step(dec, rHD);

    % Errors computation
    channelErrorStats = step(channelErrorRate, encodedMsg, rHD);
    if cnumerror >= 0
        ccHD = step(enc,guess);
        decodingErrorStats = step(decodingErrorRate, encodedMsg, ccHD);
    else
        decodingErrorStats = step(decodingErrorRate, encodedMsg, rHD);
    end

    if( sum( xor(msg,guess) ) )
        ne = ne+1;
    end
    n = n + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop(msg block) end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% logging loop time
deltaT = toc;
disp(" Simulation time: "+deltaT+" s");

% writing into results variable
results(1,1:12) = [ ...
    -1, deltaT, ...
    n, n, 0, ne, ne/n, ...
    decodingErrorStats(3), ...
    decodingErrorStats(2), decodingErrorStats(1), ...
    channelErrorStats(2), channelErrorStats(1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saving results
saveResults(output_method,results,header,["INPUT_HD_BM_BCH",method]);