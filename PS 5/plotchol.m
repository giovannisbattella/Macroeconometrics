function plotchol(cholesky, varnames, shockname, cumulate, shock, upper, lower, prc)

% Function to plot the IRFs of an SVAR identified using the Cholesky method.
% If upper and lower are provided, then also the confidence bands are
% plotted.


% Inputs:   cholesky = (N x N x horizon) array of wold impulse responses 
%           varnames = Nx1 vector of variable names
%           shock = integer indicating the shock of interest
%           upper = (N x N x horizon) array of upper bounds for IRF
%           lower = (N x N x horizon) array of lower bounds for IRF


% Returns:  plotted Wold IRF


[N,~,horizon] = size(cholesky);

irfs_plot = cholesky(:,shock,:);
irfs_plot(cumulate,:) = cumsum(irfs_plot(cumulate,:),2);

if nargin > 5
    upper_plot = upper(:,shock,:); % Assuming that these have been cumulated
    lower_plot = lower(:,shock,:);
end

nrow = ceil(sqrt(N));
ncol = ceil(N / nrow);

figure;

for i=1:N
    
    subplot(nrow, ncol, i)
    
    if nargin > 5
    
    fill([0:horizon-1 fliplr(0:horizon-1)]' ,[upper_plot(i, :)'; fliplr(lower_plot(i, :))'],...
    [0 0.4470 0.7410],'EdgeColor','None'); hold on;
    plot(0:horizon-1,(irfs_plot(i, :)),'-','LineWidth',1.5,'Color','k'); hold on;
    line(get(gca,'Xlim'),[0 0],'Color',[1 0 0],'LineStyle','--','LineWidth',1); hold off;
    
    else
        
    plot(0:horizon-1,(irfs_plot(i, :)),'-','LineWidth',1.5,'Color','k'); hold on;
    line(get(gca,'Xlim'),[0 0],'Color',[1 0 0],'LineStyle','--','LineWidth',1); hold off;
        
    end

    ylabel(varnames{i}, 'FontSize', 16);   
    title(shockname, 'FontSize', 16);  
    xlim([0 horizon-1]);
    axis tight
    set(gca,'FontSize',16)
    
end

if nargin > 5
    legend({strcat(num2str(prc), '% confidence bands'),'Point IRF'},'FontSize',16,'Orientation','Horizontal')
else
    legend('Cholesky IRF','FontSize',16,'Orientation','Horizontal')

end