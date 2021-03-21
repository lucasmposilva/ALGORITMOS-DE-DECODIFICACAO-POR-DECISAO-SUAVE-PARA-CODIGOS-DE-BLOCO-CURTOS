function FER_BER_HD_THEO(EbNoMin,EbNoMax,step,n,k,t,dmin,name)
%% Coding rate
%EbNo =; n =; k =; t =; dmin;
rc = k/n ;
results = double(zeros( (EbNoMax-EbNoMin)/step, 5));
%% Main loop
for EbNo = EbNoMin:step:EbNoMax
%%
    % dB --> linear
    EbNoLin = 10^(EbNo/10);
    
    % AWGN chanel error rate
    p       = 0.5*erfc(sqrt(EbNoLin));
    prc     = 0.5*erfc(sqrt(EbNoLin*rc));
    
    % FER
    FER = double(0);
    for j = t+1:1:n
        FER = FER + ( factorial(n)/(factorial(n-j)*factorial(j)) )*(prc^j)*((1-prc)^(n-j));
    end
    
    %BER
    BER = (dmin/n)*FER;
%%
results( round(1 + (EbNo-EbNoMin)/step), : ) = [ ...
    EbNo, FER, BER, p, prc ...
    ];
end % main loop
%%
filename = sprintf('result_BER_FER_%s_%d_%d',name,n,k);
save(filename,'results');