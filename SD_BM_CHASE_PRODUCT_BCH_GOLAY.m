function SD_BM_CHASE_PRODUCT_BCH_GOLAY(n2,k2,radius2,iterationWeight,wethr,vEbNo,output_method)

Nindexes    = size(vEbNo,2);
nIterations = size(iterationWeight,2);

% parameters
n1=23;k1=12;dmin1=7;

K = k1*k2;
N = n1*n2;

% words 1
load('words_golay.mat');
map1b0  = (words==0);
map1b1  = words;
words1  = 1 - 2*words;

% loading g, h and s matrices
load('env_GOLAY.mat')

% words 2
% They are going to be defined later by Chase

clear words;

% allocating memory
msg         = zeros(k1,k2);
msgBlock    = zeros(n1,n2);

% generator polynomials
%[gpol1, t1] = bchgenpoly(n1,k1,[],'double');
[gpol2, t2] = bchgenpoly(n2,k2,[],'double');

% encoders
%enc1 = comm.BCHEncoder(n1,k1,gpol1);
enc2 = comm.BCHEncoder(n2,k2,gpol2);

% decoders
%dec1 = comm.BCHDecoder(n1,k1,gpol1);
dec2 = comm.BCHDecoder(n2,k2,gpol2);

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
        % Tx
        % line enc
        for line = 1:1:k1
            data              = randi([0 1], k2, 1);
            msgBlock(line,:)  = step(enc2, data)';
        end
        % column enc
        for column = 1:1:n2
            msgBlock(:,column)  = encode(msgBlock(1:k1,column),n1,k1,'linear/binary',g);
        end
        
        % mod
        modSignal = 1 - 2*msgBlock;
        
        % channel
        r = step(chan, modSignal);
        
        % Rx
        rHD = (r<0);
        r   = real(r);
        
        % DEC
               
                     
        % LLR
        llrM    = r;
        for iteration = 1:nIterations
            %LINES ITERATIONS
            % lines correlation
            r2              = llrM;
            % Candidates 2
            [cand2, map2b0Multi, map2b1Multi] = prodCodeChaseCandidates(k1,n2,k2,gpol2,radius2,r2);
            
            %
            for i = 1:n2
                for j = 1:k1
                    %
                    correlation     = r2(j,:) * cand2(:,:,j);
                    map2b0          = map2b0Multi(:,:,j);
                    map2b1          = map2b1Multi(:,:,j);
                    %
                    if( sum(map2b0(i,:))~=0 && sum(map2b1(i,:))~=0 )
                        %
                        [m0, key0]  = max(map2b0(i,:).*correlation);
                        [m1, key1]  = max(map2b1(i,:).*correlation);
                        %
                        s = sum( r2(j,:).*(cand2(:,key0,j)' - cand2(:,key1,j)')./2 );
                        %
                        c = r2(j,i).*(cand2(i,key0,j)' - cand2(i,key1,j)')./2;  
                        %
                        r2(j,i) = r2(j,i) + iterationWeight(iteration).*(s - c);
                        %
                    else
                        if( sum(map2b0(i,:))==0 )
                            r2(j,i) = r2(j,i) - iterationWeight(iteration);
                        else
                            r2(j,i) = r2(j,i) + iterationWeight(iteration);
                        end
                    end
                end
            end
            % UPDATING
            llrM = r2;
            
            %COLUMNS ITERATIONS
            % columns correlation
            r1             = llrM';
            correlation    = r1 * words1;
            %
            for i = 1:n1
                for j = 1:n2
                    %
                    [m0, key0]  = max(map1b0(i,:).*correlation(j,:));
                    [m1, key1]  = max(map1b1(i,:).*correlation(j,:));
                    %
                    s = sum( r1(j,:).*( words1(:,key0)' - words1(:,key1)' )./2 );
                    %
                    c = r1(j,i).*( words1(i,key0)' - words1(i,key1)' )./2;  
                    %
                    r1(j,i) = r1(j,i) + iterationWeight(iteration).*(s - c);
                    %
                end
            end
            % UPDATING
            llrM = r1';
        end
        %
        choice = prodCodeDec(n1,n2,words1,llrM');
        % Mapping +1 and -1 to 0 and 1
        choice = ((1 - choice)./2)';
        % Errors computation
        for x = 1:n1
            channelErrorStats   = step(channelErrorRate, msgBlock(x,:)', double(rHD(x,:)'));
            decodingErrorStats  = step(decodingErrorRate, msgBlock(x,:)', choice(x,:)');
        end
        if( sum(sum(xor(msgBlock',choice'))) )
            ne = ne + 1;
        end
        n = n + 1;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % logging loop time
    deltaT = toc;
    disp("EbNo: "+EbNo+"dB"+" Simulation time: "+deltaT+" s")
    
    % writing into results variable
    results(INDEX,1:12) = [ ...
        vEbNo(INDEX), deltaT, ...
        n, 0, n, ...
        ne, ne/n, ...
        decodingErrorStats(3), ...
        decodingErrorStats(2), decodingErrorStats(1), ...
        channelErrorStats(2), channelErrorStats(1), ...
        ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ROUND end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saving results
saveResults(output_method,results,header,["SD_BM_CHASE_PRODUCT_BCH",mat2str(vEbNo)]);