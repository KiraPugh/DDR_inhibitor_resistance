%% 1. IMPORT THE DATA

% Define parameters
dose_increments = 0:10:100;  % Dose increments doses of drugs = 0,10,20,...,100 
no_circles = [1, 2, 4, 8, 16, 32, 64, 128, 256, 482];  % Different numbers of clusters (mono to single-cell)

% Initialize a cell array to hold all data
dataCells = cell(length(no_circles), length(dose_increments));

% Loop through each NoCircles value
for n = 1:length(no_circles)
    % Loop through each dose increment
    for d = 1:length(dose_increments)
        % Generate the filename
        filename = sprintf('oct_B_mu41000_sigma8200_DoseC%d_DoseO%d_NoCircles%d_DRfrac30_DRtype3', ...
            dose_increments(d), dose_increments(d), no_circles(n));
        
        % Read the table and assign it to the cell array
        dataCells{n, d} = readtable(filename);
    end
end

%% 2. SORT OUT THE INDIVIDUAL RUNS OF THE DATA

% Define the number of runs and the number of rows per run
numRuns = 100; %number of times we run the in silico experiments
rowsPerRun = 311;

% Initialize cell array to store processed data for each heatmap point
cells_totcells = cell(size(dataCells, 1), size(dataCells, 2), numRuns);
cells_DRfrac = cell(size(dataCells, 1), size(dataCells, 2), numRuns);


for n = 1:size(dataCells, 1)
    for d = 1:size(dataCells, 2)
        for run = 1:numRuns
            startRow = (run - 1) * rowsPerRun + 1;
            endRow = run * rowsPerRun;
            
            % Extract the data for the current run
            cells_totcells{n, d, run} = dataCells{n, d}{startRow:endRow, 5};  % total number of cells
            cells_DRfrac{n, d, run} = dataCells{n, d}{startRow:endRow, 15} ./ dataCells{n, d}{startRow:endRow, 5};  % fraction of drug resistant cells
            
            
        end
    end
end

% Calculate the mean for each heatmap point across runs
means_totcells = cell(size(dataCells));
means_DRfrac = cell(size(dataCells));

for n = 1:size(dataCells, 1)
    for d = 1:size(dataCells, 2)
        % Concatenate data for each cell and calculate the mean
        stack_totcells = cat(3, cells_totcells{n, d, :});
        stack_DRfrac = cat(3, cells_DRfrac{n, d, :});
        means_totcells{n, d} = mean(stack_totcells, 3);
        means_DRfrac{n, d} = mean(stack_DRfrac, 3);
    end
end

%% 3. PLOT GRAPH

% set up a matrix of zeros for the heatmaps
data_for_heatmap_totcells = zeros(size(dataCells, 1), size(dataCells, 2));
data_for_heatmap_DRfrac = zeros(size(dataCells, 1), size(dataCells, 2));

% Fill the data_for_heatmap matrix with the last value from the mean of each cell
for n = 1:size(dataCells, 1)
    for d = 1:size(dataCells, 2)
        data_for_heatmap_totcells(n, d) = means_totcells{n, d}(end);
        data_for_heatmap_DRfrac(n, d) = means_DRfrac{n, d}(end);
    end
end


%get matrix in correct order
data_for_heatmap = flipud(fliplr(data_for_heatmap'));
data_for_heatmap2 = flipud(fliplr(data_for_heatmap2'));


figure(1)

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

