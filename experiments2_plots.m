%% 1. IMPORT THE DATA

% Define parameters
dose_increments = 0:10:100;  % Dose increments
no_circles = [1, 2, 4, 8, 16, 32, 64, 128, 256, 482];  % Different NoCircles values

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
numRuns = 100;
rowsPerRun = 311;

% Initialize cell array to store processed data for each heatmap point
cells = cell(size(dataCells, 1), size(dataCells, 2), numRuns);
cells2 = cell(size(dataCells, 1), size(dataCells, 2), numRuns);


% Loop through each NoCircles value
for n = 1:size(dataCells, 1)
    % Loop through each dose increment
    for d = 1:size(dataCells, 2)
        % Loop through the runs
        for run = 1:numRuns
            startRow = (run - 1) * rowsPerRun + 1;
            endRow = run * rowsPerRun;
            
            % Extract the data for the current run (e.g., 5th column)
            cells{n, d, run} = dataCells{n, d}{startRow:endRow, 5};  % total number of cells
            cells2{n, d, run} = dataCells{n, d}{startRow:endRow, 15} ./ dataCells{n, d}{startRow:endRow, 5};  % fraction of drug resistant cells


        end
    end
end

% Calculate the mean for each heatmap point across runs
means = cell(size(dataCells));
means2 = cell(size(dataCells));

for n = 1:size(dataCells, 1)
    for d = 1:size(dataCells, 2)
        % Concatenate data for each cell and calculate the mean
        stack = cat(3, cells{n, d, :});
        stack2 = cat(3, cells2{n, d, :});
        means{n, d} = mean(stack, 3);
        means2{n, d} = mean(stack2, 3);
    end
end

%% 3. PLOT GRAPH

% Preallocate the data_for_heatmap matrix
data_for_heatmap = zeros(size(dataCells, 1), size(dataCells, 2));
data_for_heatmap2 = zeros(size(dataCells, 1), size(dataCells, 2));

% Fill the data_for_heatmap matrix with the last value from the mean of each cell
for n = 1:size(dataCells, 1)
    for d = 1:size(dataCells, 2)
        data_for_heatmap(n, d) = means{n, d}(end);
        data_for_heatmap2(n, d) = means2{n, d}(end);  
    end
end

data_for_heatmap = flipud(fliplr(data_for_heatmap'));
data_for_heatmap2 = flipud(fliplr(data_for_heatmap2'));


% Create a figure with subplots
figure(1)



% Subplot 1 for data_for_heatmap
subplot(1, 2, 1);  % 1 row, 2 columns, 1st subplot
imagesc(data_for_heatmap);
colormap('copper');
colorbar;
xLabels = {'482','256','128','64','32','16','8','4','2','1'}; % X-axis labels
yLabels = {'1','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1','0'}; % Y-axis labels
set(gca, 'XTick', 1:size(data_for_heatmap, 2), 'XTickLabel', xLabels);
set(gca, 'YTick', 1:size(data_for_heatmap, 1), 'YTickLabel', yLabels);
xlabel('Number of seeded clusters');
ylabel('Dose of drugs 1 and 2 (\muM)');
title('Total Cell Count');

% Subplot 2 for data_for_heatmap2
subplot(1, 2, 2);  % 1 row, 2 columns, 2nd subplot
imagesc(data_for_heatmap2);
colormap('parula');
colorbar;
xLabels = {'482','256','128','64','32','16','8','4','2','1'}; % X-axis labels
yLabels = {'1','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1','0'}; % Y-axis labels
set(gca, 'XTick', 1:size(data_for_heatmap2, 2), 'XTickLabel', xLabels);
set(gca, 'YTick', 1:size(data_for_heatmap2, 1), 'YTickLabel', yLabels);
xlabel('Number of seeded clusters');
ylabel('Dose of drugs 1 and 2 (\muM)');
title('Fraction of Resistant Cells');

