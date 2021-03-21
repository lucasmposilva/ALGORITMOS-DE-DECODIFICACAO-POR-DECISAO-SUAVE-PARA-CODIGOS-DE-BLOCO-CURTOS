function resultsToPlot(file,xArr,yArr,titles,lineTypes,legends)
%

% loading results variable from file
load(file,'results');

% calling a new figure
figure(1); clf;

% drawing curves
for key = 1:size(xArr,2)
    cmd = "semilogy(results(%s),results(%s),'%s','LineWidth',1,'MarkerSize',4)";
    cmd = sprintf( cmd, xArr(key), yArr(key), lineTypes(key) );
    eval(cmd);
    hold on;
end

hold off;

% titles
title( { strrep( titles(1), "_", " "), ''} );
xlabel(titles(2));
ylabel(titles(3));
legend(legends,'Location','northeast');

grid on;
set(gca,'YMinorTick','on');
pbaspect([3 4 1]);

% writing a file
print(1,titles(1),'-dpng','-r400');