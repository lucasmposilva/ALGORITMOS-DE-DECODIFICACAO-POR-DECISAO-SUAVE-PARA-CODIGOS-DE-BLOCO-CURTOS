function [choice] = prodCodeDec(n1,n2,words,llrM)
    choice      = zeros(n2,n1); 
    correlation = llrM * words;
    [maxV, maxK]= max(correlation');
    
    choice(1:size(maxK,2),:) = words(:,maxK)';
    
    

