function [status, choice1] = prodCodeDecMatch(k1,n1,k2,n2,words1,words2,llrM)
    %
    status = 0;
    
    % 1
    choice1      = zeros(n2,n1); 
    correlation = llrM * words1;
    [maxV, maxK]= max(correlation');
    
    choice1(1:size(maxK,2),:) = words1(:,maxK)';

    % 2
    choice2      = zeros(n1,n2); 
    correlation = llrM' * words2;
    [maxV, maxK]= max(correlation');
    
    choice2(1:size(maxK,2),:) = words2(:,maxK)';
    
    % compare
    if( sum(sum(xor(choice1(1:n2,1:k1)',choice2(1:k1,1:n2)))) == 0 )
        status = 1;
    end