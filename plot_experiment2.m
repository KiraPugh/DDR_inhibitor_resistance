%load the C++ output
load('output_exp2.mat');


%% PLOT HEATMAPS

data_for_heatmap_totcells = zeros(size(dataCells, 1), size(dataCells, 2));
data_for_heatmap_DRfrac = zeros(size(dataCells, 1), size(dataCells, 2));

% Fill the data_for_heatmap matrix with the last value from the mean of each cell
for n = 1:size(dataCells, 1)
    for d = 1:size(dataCells, 2)
        data_for_heatmap_totcells(n, d) = means_totcells{n, d}(end);
        data_for_heatmap_DRfrac(n, d) = means_DRfrac{n, d}(end);
    end
end

figure;

% Subplot 1 for data_for_heatmap_totcells
subplot(1, 2, 1);
imagesc(data_for_heatmap_totcells);
colormap('copper');
colorbar;
xLabels = {'482','256','128','64','32','16','8','4','2','1'};
yLabels = {'1','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1','0'};
set(gca, 'XTick', 1:size(data_for_heatmap_totcells, 2), 'XTickLabel', xLabels);
set(gca, 'YTick', 1:size(data_for_heatmap_totcells, 1), 'YTickLabel', yLabels);
xlabel('Number of seeded clusters');
ylabel('Dose of drugs 1 and 2 (\muM)');
title('Total Cell Count');

% Subplot 2 for data_for_heatmap_DRfrac
subplot(1, 2, 2);
imagesc(data_for_heatmap_DRfrac);
colormap('parula');
colorbar;
xLabels = {'482','256','128','64','32','16','8','4','2','1'};
yLabels = {'1','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1','0'};
set(gca, 'XTick', 1:size(data_for_heatmap_DRfrac, 2), 'XTickLabel', xLabels);
set(gca, 'YTick', 1:size(data_for_heatmap_DRfrac, 1), 'YTickLabel', yLabels);
xlabel('Number of seeded clusters');
ylabel('Dose of drugs 1 and 2 (\muM)');
title('Fraction of Drug Resistant Cells');


%% PLOT DYNAMIC TOTAL NUMBER OF CELLS

% Define the time vector
time = [0:1:310]';

% Number of rows and columns for the subplot grid
numRows = size(dataCells, 2); % number of rows of data
numCols = size(dataCells, 1); % number of columns of data

xLabels = {'482','256','128','64','32','16','8','4','2','1'};
yLabels = {'1','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1','0'};

figure;

for n = 1:size(dataCells, 1)
    for d = 1:size(dataCells, 2)
        % Calculate the index for the subplot
        subplotIndex = (n-1) * size(dataCells, 2) + d;
        
        % Create a subplot in the current index
        subplot(size(dataCells, 1), size(dataCells, 2), subplotIndex);
        
        % Extract the mean values for the current cell
        mean_values = means_totcells{n, d};
        
        % Plot the mean values over time
        plot(time, mean_values, 'LineWidth', 2);
        
        % Set title, labels, and limits
        title(sprintf('DR Fraction (%s,%s)', xLabels{d}, yLabels{n}));
        xlabel('Time');
        ylabel('Cell Count');
        xlim([0 310]);
    end
end
sgtitle('Dynamic Cell Count');

%% PLOT DYNAMIC FRACTION OF DRUG RESISTANT CELLS

% Define the time vector
time = [0:1:310]';

figure;

for n = 1:size(dataCells, 1)
    for d = 1:size(dataCells, 2)
        % Calculate the index for the subplot
        subplotIndex = (n-1) * size(dataCells, 2) + d;
        
        % Create a subplot in the current index
        subplot(size(dataCells, 1), size(dataCells, 2), subplotIndex);
        
        % Extract the mean values for the current cell
        mean_values = means_DRfrac{n, d};
        
        % Plot the mean values over time
        plot(time, mean_values, 'LineWidth', 2);
        
        % Set the title, labels, and limits
        title(sprintf('DR Fraction (%s,%s)', xLabels{d}, yLabels{n}));
        xlabel('Time');
        ylabel('DR Frac');
        xlim([0 310]);
    end
end
sgtitle('Dynamic DR Fraction for (No Clusters, Drug Dose)');
