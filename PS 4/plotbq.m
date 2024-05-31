function plotLR(irf, varnames, shocknames, cumulate, shocks, upper, lower, prc)

% Function to plot the IRFs of an SVAR identified using the LR method.
% If upper and lower are provided, then also the confidence bands are
% plotted.


% Inputs:   irf = (N x N x horizon) array of LR impulse responses 
%           varnames = Nx1 vector of variable names
%           shocks = integer indicating the shock(s) of interest
%           upper = (N x N x horizon) array of upper bounds for IRF
%           lower = (N x N x horizon) array of lower bounds for IRF
%           prc     = integer confidence bands width

%%
% Returns:  plotted LR IRF

[N,~,horizon] = size(irf);

irfs_plot = irf(:,shocks,:);

for s=shocks
    irfs_plot(cumulate,s,:) = cumsum(irfs_plot(cumulate,s,:),3);
end

if nargin > 5
    upper_plot = upper(:,shocks,:); % Assuming that these have been cumulated
    lower_plot = lower(:,shocks,:);
end


nrow = size(varnames,2);
ncol = size(shocknames,2);

counter  = 0;

figure;
for i=1:nrow
    for j=1:ncol
    
        counter = counter + 1;
        subplot(nrow, ncol, counter)

        if nargin > 5

        fill([0:horizon-1 fliplr(0:horizon-1)]' ,[upper_plot(i,j,:)'; fliplr(lower_plot(i,j,:))'],...
        [0 0.4470 0.7410],'EdgeColor','None'); hold on;
        plot(0:horizon-1,(irfs_plot(i,j,:)),'-','LineWidth',1.5,'Color','k'); hold on;
        line(get(gca,'Xlim'),[0 0],'Color',[1 0 0],'LineStyle','--','LineWidth',1); hold off;

        else

        plot(0:horizon-1,(irfs_plot(i,j,:)),'-','LineWidth',1.5,'Color','k'); hold on;
        line(get(gca,'Xlim'),[0 0],'Color',[1 0 0],'LineStyle','--','LineWidth',1); hold off;

        end

        ylabel(varnames{i}, 'FontSize', 16);   
        title(shocknames{j}, 'FontSize', 16);  
        xlim([0 horizon-1]);
        axis tight
        set(gca,'FontSize',16)
    end
end

if nargin > 5
    legend({strcat(num2str(prc), '% confidence bands'),'Point IRF'},'FontSize',16,'Orientation','Horizontal')
else
    legend('Point IRF','FontSize',16,'Orientation','Horizontal')

end