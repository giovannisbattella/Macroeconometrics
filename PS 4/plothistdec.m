function plothistdec(histdec,ystar,shocks,series,dates,shocknames,varnames)

% Function to plot the historical decomposition of the series of choice to
% the shocks of choice

% Inputs:   histdec     = T-p x N x N matrix of contributions of shocks to series
%           ystar       = T-p x N demeaned series
%           shocks       = position of shocks of interest
%           series      = vector of chosen series
%           shocknames   = names of chosen shocks
%           varnames    = names of variables in chosen series

% Returns: plots of the chosen series to the chosen shock

nplots = size(series,2);
figure;
for i=1:nplots
    
    subplot(nplots,1,i)
    bar(dates,squeeze(histdec(:,:,i)), 'stacked'); hold on;
    plot(dates,ystar(:,i),'LineWidth',1.25, 'Color', 'k');  hold off;
    recessionplot
    axis tight
    set(gca,'FontSize',15)
    title(varnames{i},'Interpreter','Latex')
    xlim([min(dates), max(dates)])
       
end
legend([shocknames, 'Series'], 'Orientation','horizontal','Interpreter','Latex')

end