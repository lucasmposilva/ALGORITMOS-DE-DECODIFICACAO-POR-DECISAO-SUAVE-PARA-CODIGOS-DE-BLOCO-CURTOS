function saveResults(mode,results,header,log)
%
%function saveResults(mode,results,header,log)
%
%Save the variables 'results', 'header' and 'log' into a file called 'results.mat'.
%
%   mode:
%       'w'   - write mode
%       'a'   - append mode
%
%   results:
%       a variable that contains the simulation results.
%
%   header:
%       a variable that describes each column of the variable 'results'.
%       e.g.: ["EbNo","Bits","Errors","BER"]
%
%   log:
%       a variable that describes the simulation, i.e. the order
%       that the scripts were called and for what values of EbNo they were called.
%       e.g.: ["scriptA","[1,2,3,4]";"scriptB","[0,1]"]
%       note: You have to use double quotes.
%

% write mode
if(mode=='w')
    save('results.mat','results','header','log');
% append mode
elseif(mode=='a')
    % copying
    r = results; h = header; l = log;
    % loading
    load('results.mat');
    % appending
    results = [results; r];
    if (sum(h == header) ~= size(header,2))
        disp('Error: unable to append results with different headers');
        return;
    end
    log     = [log; l];
    % saving
    save('results.mat','results','header','log');
% wrong mode
else
    disp('Wrong call to function saveResults...');
    help saveResults;
end
%
end