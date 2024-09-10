%load the output from C++
load('output_exp3.mat');

%% PLOT HEATMAPS

data_for_heatmap_totcells_singlecellclusters = zeros(size(dataCells, 1), size(dataCells, 2));
data_for_heatmap_totcells_multicellclusters = zeros(size(dataCells, 1), size(dataCells, 2));
data_for_heatmap_totcells_monoclusters = zeros(size(dataCells, 1), size(dataCells, 2));
data_for_heatmap_DRfrac_singlecellclusters = zeros(size(dataCells, 1), size(dataCells, 2));
data_for_heatmap_DRfrac_multicellclusters = zeros(size(dataCells, 1), size(dataCells, 2));
data_for_heatmap_DRfrac_monoclusters = zeros(size(dataCells, 1), size(dataCells, 2));


% Fill the data_for_heatmap matrix with the last value from the mean of each cell
for n = 1:size(dataCells, 1)
    for d = 1:size(dataCells, 2)
        data_for_heatmap_totcells_singlecellclusters(n, d) = means_totcells_singlecellclusters{n, d}(end);
        data_for_heatmap_totcells_multicellclusters(n, d) = means_totcells_multicellclusters{n, d}(end);
        data_for_heatmap_totcells_monoclusters(n, d) = means_totcells_monoclusters{n, d}(end);
        data_for_heatmap_DRfrac_singlecellclusters(n, d) = means_DRfrac_singlecellclusters{n, d}(end);
        data_for_heatmap_DRfrac_multicellclusters(n, d) = means_DRfrac_multicellclusters{n, d}(end);
        data_for_heatmap_DRfrac_monoclusters(n, d) = means_DRfrac_monoclusters{n, d}(end);
    end
end

% Subplot 1: Single-cell Total Cells
subplot(2, 3, 1);
imagesc(data_for_heatmap_totcells_singlecellclusters);
colormap('copper');
colorbar;
set(gca, 'XTick', 1:size(data_for_heatmap_totcells_singlecellclusters, 2), ...
    'XTickLabel', {'41','49.2','57.4','65.6','73.8','82'}, ...
    'YTick', 1:size(data_for_heatmap_totcells_singlecellclusters, 1), ...
    'YTickLabel', {'1','0.8','0.6','0.4','0.2','0'});
xlabel('\tau_R');
ylabel('Dose of drugs 1 and 2 (\muM)');
title('Single-cell clusters: Total cell count');

% Subplot 2: Multi-cell Total Cells
subplot(2, 3, 2);
imagesc(data_for_heatmap_totcells_multicellclusters);
colormap('copper');
colorbar;
set(gca, 'XTick', 1:size(data_for_heatmap_totcells_multicellclusters, 2), ...
    'XTickLabel', {'41','49.2','57.4','65.6','73.8','82'}, ...
    'YTick', 1:size(data_for_heatmap_totcells_multicellclusters, 1), ...
    'YTickLabel', {'1','0.8','0.6','0.4','0.2','0'});
xlabel('\tau_R');
ylabel('Dose of drugs 1 and 2 (\muM)');
title('Multi-cell clusters: Total cell count');

% Subplot 3: Mono-clusters Total Cells
subplot(2, 3, 3);
imagesc(data_for_heatmap_totcells_monoclusters);
colormap('copper');
colorbar;
set(gca, 'XTick', 1:size(data_for_heatmap_totcells_monoclusters, 2), ...
    'XTickLabel', {'41','49.2','57.4','65.6','73.8','82'}, ...
    'YTick', 1:size(data_for_heatmap_totcells_monoclusters, 1), ...
    'YTickLabel', {'1','0.8','0.6','0.4','0.2','0'});
xlabel('\tau_R');
ylabel('Dose of drugs 1 and 2 (\muM)');
title('Monoclusters: Total cell count');

% Subplot 4: Single-cell Fraction DR
subplot(2, 3, 4);
imagesc(data_for_heatmap_DRfrac_singlecellclusters);
colormap('copper');
colorbar;
set(gca, 'XTick', 1:size(data_for_heatmap_DRfrac_singlecellclusters, 2), ...
    'XTickLabel', {'41','49.2','57.4','65.6','73.8','82'}, ...
    'YTick', 1:size(data_for_heatmap_DRfrac_singlecellclusters, 1), ...
    'YTickLabel', {'1','0.8','0.6','0.4','0.2','0'});
xlabel('\tau_R');
ylabel('Dose of drugs 1 and 2 (\muM)');
title('Single-cell clusters: Fraction of drug-resistant cells');

% Subplot 5: Multi-cell Fraction DR
subplot(2, 3, 5);
imagesc(data_for_heatmap_DRfrac_multicellclusters);
colormap('copper');
colorbar;
set(gca, 'XTick', 1:size(data_for_heatmap_DRfrac_multicellclusters, 2), ...
    'XTickLabel', {'41','49.2','57.4','65.6','73.8','82'}, ...
    'YTick', 1:size(data_for_heatmap_DRfrac_multicellclusters, 1), ...
    'YTickLabel', {'1','0.8','0.6','0.4','0.2','0'});
xlabel('\tau_R');
ylabel('Dose of drugs 1 and 2 (\muM)');
title('Multi-cell clusters: Fraction of drug-resistant cells');

% Subplot 6: Mono-clusters Fraction DR
subplot(2, 3, 6);
imagesc(data_for_heatmap_DRfrac_monoclusters);
colormap('copper');
colorbar;
set(gca, 'XTick', 1:size(data_for_heatmap_DRfrac_monoclusters, 2), ...
    'XTickLabel', {'41','49.2','57.4','65.6','73.8','82'}, ...
    'YTick', 1:size(data_for_heatmap_DRfrac_monoclusters, 1), ...
    'YTickLabel', {'1','0.8','0.6','0.4','0.2','0'});
xlabel('\tau_R');
ylabel('Dose of drugs 1 and 2 (\muM)');
title('Monoclusters: Fraction of drug-resistant cells');





xLabels = {'41','49.2','57.4','65.6','73.8','82'};
yLabels = {'1','0.8','0.6','0.4','0.2','0'};


%% PLOT DYNAMIC CELL COUNT (single-cell clusters)

% Define the time vector
time = [0:1:310]';

% Number of cells in each dimension
numRows = size(dataCells, 1);
numCols = size(dataCells, 2);

figure;

for n = 1:numRows
    for d = 1:numCols
        % Calculate the index for the subplot
        subplotIndex = (n - 1) * numCols + d;
        
        % Create a subplot for the current index
        subplot(numRows, numCols, subplotIndex);
        
        % Extract the mean values for the current cell
        mean_values = means_totcells_singlecellclusters{n, d};
        
        % Plot the mean values over time
        plot(time, mean_values, 'LineWidth', 2);
        
        % Set subplot title and labels
        title(sprintf('TotCells (\\tau_R %s, dose %s)', xLabels{d}, yLabels{n}));
        xlabel('Time');
        ylabel('Cell Count');
        xlim([0 310]);
    end
end
sgtitle('Dynamic Cell Count for Single-cell Clusters');





%% PLOT DYNAMIC CELL COUNT (multi-cell clusters)

% Define the time vector
time = [0:1:310]';

% Number of cells in each dimension
numRows = size(dataCells, 1);
numCols = size(dataCells, 2);

figure;

for n = 1:numRows
    for d = 1:numCols
        % Calculate the index for the subplot
        subplotIndex = (n - 1) * numCols + d;
        
        % Create a subplot for the current index
        subplot(numRows, numCols, subplotIndex);
        
        % Extract the mean values for the current cell
        mean_values = means_totcells_multicellclusters{n, d};
        
        % Plot the mean values over time
        plot(time, mean_values, 'LineWidth', 2);
        
        title(sprintf('TotCells (\\tau_R %s, dose %s)', xLabels{d}, yLabels{n}));
        xlabel('Time');
        ylabel('Cell Count');
        xlim([0 310]);
        
        grid off;
    end
end
sgtitle('Dynamic Cell Count for Multi-cell Clusters');



%% PLOT DYNAMIC CELL COUNT (monoclusters)

% Define the time vector
time = [0:1:310]';

% Number of cells in each dimension
numRows = size(dataCells, 1);
numCols = size(dataCells, 2);

figure;

for n = 1:numRows
    for d = 1:numCols
        % Calculate the index for the subplot
        subplotIndex = (n - 1) * numCols + d;
        
        % Create a subplot for the current index
        subplot(numRows, numCols, subplotIndex);
        
        % Extract the mean values for the current cell
        mean_values = means_totcells_monoclusters{n, d};
        
        % Plot the mean values over time
        plot(time, mean_values, 'LineWidth', 2);
        
        title(sprintf('TotCells (\\tau_R %s, dose %s)', xLabels{d}, yLabels{n}));
        xlabel('Time');
        ylabel('Cell Count');
        xlim([0 310]);
    end
end
sgtitle('Dynamic Cell Count for Monoclusters');





%% PLOT DYNAMIC DR FRAC (single-cell clusters)

% Define the time vector
time = [0:1:310]';

% Number of cells in each dimension
numRows = size(dataCells, 1);
numCols = size(dataCells, 2);

figure;

for n = 1:numRows
    for d = 1:numCols
        % Calculate the index for the subplot
        subplotIndex = (n - 1) * numCols + d;
        
        % Create a subplot for the current index
        subplot(numRows, numCols, subplotIndex);
        
        % Extract the mean values for the current cell
        mean_values = means_DRfrac_singlecellclusters{n, d};
        
        % Plot the mean values over time
        plot(time, mean_values, 'LineWidth', 2);
        
        title(sprintf('TotCells (\\tau_R %s, dose %s)', xLabels{d}, yLabels{n}));
        xlabel('Time');
        ylabel('DR Fraction');
        xlim([0 310]);
    end
end
sgtitle('Dynamic DR Fraction for Single-cell Clusters');


%% PLOT DYNAMIC DR FRAC (multi-cell clusters)
% Define the time vector
time = [0:1:310]';

% Number of cells in each dimension
numRows = size(dataCells, 1); % Get the number of rows
numCols = size(dataCells, 2); % Get the number of columns

figure;

for n = 1:numRows
    for d = 1:numCols
        % Calculate the index for the subplot
        subplotIndex = (n - 1) * numCols + d;
        
        % Create a subplot for the current index
        subplot(numRows, numCols, subplotIndex);
        
        % Extract the mean values for the current cell
        mean_values = means_DRfrac_multicellclusters{n, d};
        
        % Plot the mean values over time
        plot(time, mean_values, 'LineWidth', 2);
        
        title(sprintf('TotCells (\\tau_R %s, dose %s)', xLabels{d}, yLabels{n}));
        xlabel('Time');
        ylabel('DR Fraction');
        xlim([0 310]);
    end
end
sgtitle('Dynamic DR Fraction for Multi-cell Clusters');



%% PLOT DYNAMIC DR FRAC (monoclusters)

% Define the time vector
time = [0:1:310]';

% Number of cells in each dimension
numRows = size(dataCells, 1);
numCols = size(dataCells, 2);

figure;

for n = 1:numRows
    for d = 1:numCols
        % Calculate the index for the subplot
        subplotIndex = (n - 1) * numCols + d;
        
        % Create a subplot for the current index
        subplot(numRows, numCols, subplotIndex);
        
        % Extract the mean values for the current cell
        mean_values = means_DRfrac_monoclusters{n, d};
        
        % Plot the mean values over time
        plot(time, mean_values, 'LineWidth', 2);
        
        % Set subplot title and labels
        title(sprintf('DR Fraction (\\tau_R %s, dose %s)', xLabels{d}, yLabels{n}));
        xlabel('Time');
        ylabel('DR Fraction');
        xlim([0 310]);
    end
end
sgtitle('Dynamic DR Fraction for Monoclusters');
