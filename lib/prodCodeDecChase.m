function [choice] = prodCodeDecChase(n1,n2,words,llrM)
    choice      = zeros(n2,n1);
    
    for index = 1:size(llrM,1)
        correlation = llrM(index,:) * words(:,:,index);
        
        [maxV, maxK]= max(correlation');
    
        choice(index,:) = words(:,maxK,index)';
    end
    

