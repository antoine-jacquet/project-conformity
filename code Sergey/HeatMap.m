clear
close all

load('ConformitySimulations.mat');
Data = ConformitySimulations;

[~, ind] = unique(ConformitySimulations(:, 1:10), 'rows');

Data = ConformitySimulations(ind, :);

ntot = Data(:,1);
ind = find(ntot==1000);
Data = Data(ind, :);

%% COLOR PLOT FOR SUCCESS RATE %%
%for reference on variable order: Success_rate = [ntot N T T0 vs_setprefs group_invasion q0 d beta r Successes/Runs]

ntot = Data(:,1);
N = Data(:,2);
T = Data(:,3);
T0 = Data(:,4);
vs_setprefs = Data(:,5);
group_invasion = Data(:,6);
q0 = Data(:,7);
d = Data(:,8);
beta = Data(:,9);
r = Data(:,10);
success_rate = Data(:,11);

cdivs = 10;
edges = 0:1/cdivs:1;
[Count, edges, bin] = histcounts(success_rate,edges);


cmap = ones(cdivs, 3);  % RGB colors
    lowc = .1;
    highc = .9;
    step = (highc - lowc)/(cdivs-1);
    % Red gradient
    %cmap(:,2) = fliplr(lowc:step:highc);
    %cmap(:,3) = fliplr(lowc:step:highc);
    % Blue to red gradient
    %cmap(:,1) = lowc:step:highc;
    %cmap(:,2) = 0;
    %cmap(:,3) = fliplr(lowc:step:highc);
    % Yellow to green gradient
    %cmap(:,1) = fliplr(lowc:step:highc);
    %cmap(:,2) = 1;
    %cmap(:,3) = 0;
    % Red to Yellow to Green gradient
    %cmap(:,1) = [0:2/cdivs:.8 ones(1,5)];
    %cmap(:,2) = [ones(1,5) 0.8:-2/cdivs:0];
    %cmap(:,3) = 0;
    % Red to Green gradient
    %cmap(:,1) = fliplr(lowc:step:highc);
    %cmap(:,2) = lowc:step:highc;
    %cmap(:,3) = 0;
    % Red to Yellow to Green gradient
    cmap(:,1) = 0;
    cmap(:,2) = fliplr([ones(1,5) 0.8:-2/cdivs:0]);
    cmap(:,3) = fliplr([0:2/cdivs:.8 ones(1,5)]);

    R = unique(r);

for ii=1:length(R)
    
    figure;
    hold on;
    
    rTemp = R(ii);
    
    ind = find(r == rTemp);
 
    x = d(ind);
    y = beta(ind);
    
    for i=1:cdivs
       idx = bin(ind)==i;
       plot(x(idx), y(idx), '.', 'MarkerSize', 5, 'Color', cmap(i,:));
    end
    
    title({'Successful invasion of conformity (in % of all simulations), q0 = .2'; ['Recombination rate: r = ', num2str(rTemp)]})
    %subtitle(sprintf('Recombination rate: r = %0.2f', rTemp))
    ylabel('Conformity strength \beta')
    ylim([min(y)-.1 3+.1])
    xlabel('Dispersal rate d')
    xlim([-0.005 .105])

    colormap(cmap)
    caxis([0 1])
    colorbar
    
    if ismac
        png_file = sprintf('/Users/Antoine/Documents/GitHub/project-conformity/code Sergey/Figs/heatmap_r%1.2f.png',rTemp);
        pdf_file = sprintf('/Users/Antoine/Documents/GitHub/project-conformity/code Sergey/Figs/heatmap_r%1.2f.pdf',rTemp);
        eps_file = sprintf('/Users/Antoine/Documents/GitHub/project-conformity/code Sergey/Figs/heatmap_r%1.2f.eps',rTemp);
    elseif ispc
        png_file = sprintf('Z:/Conformity/MATLAB/Figs/heatmap_r%1.2f.png',rTemp);
    end
	print(png_file,'-dpng');
    print('-depsc2',eps_file);
    print('-dpdf',pdf_file,'-bestfit');

end



