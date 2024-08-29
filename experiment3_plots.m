%% 1. IMPORT THE DATA

tot_heatmappoints = 36;

% Define parameters
dose_increments = 0:20:100;  % Dose increments doses of 0,20,40,60,80,100
doubling_time_resistant = [41000, 49200, 57399, 65600, 73800, 82000];  % Doubling times of the drug resistant cells
doubling_time_resistant_std = [8200, 9840, 11480, 13120, 14760, 16400];  % Standard deviations of the drug resistant cells

% Initialize cell arrays to hold all data
dataCells = cell(3, length(doubling_time_resistant), length(dose_increments));
dataTypes = {'monoclusters', 'multicell', 'singlecell'};
noCircles = [1, 32, 419];

for t = 1:3
    for m = 1:length(doubling_time_resistant_std)
        for d = 1:length(dose_increments)
            % Generate the filename
            filename = sprintf('oct_B_mu_ds41000_sigma_ds8200_mu_dr%d_sigma_dr%d_DoseC%d_DoseO%d_NoCircles%d_DRfrac30_DRtype3', ...
                doubling_time_resistant(m), doubling_time_resistant_std(m), dose_increments(d), dose_increments(d), noCircles(t));
            
            % Read data and store in cell array
            dataCells{t, m, d} = readtable(filename);
        end
    end
end

%% 2. SORT OUT THE INDIVIDUAL RUNS OF THE DATA

numRuns = 100; % number of times we run the in silico experiments
rowsPerRun = 311; % number of time points we save (every hour for 0-310 hours)

% Initialize cell arrays for processed data
cells_totcells = cell(3, length(doubling_time_resistant), length(dose_increments), numRuns);
cells_fracDR = cell(3, length(doubling_time_resistant), length(dose_increments), numRuns);

for t = 1:3
    for n = 1:length(doubling_time_resistant_std)
        for d = 1:length(dose_increments)
            for run = 1:numRuns
                startRow = (run - 1) * rowsPerRun + 1;
                endRow = run * rowsPerRun;
                
                % Extract the data for the current run
                cells_totcells{t, n, d, run} = dataCells{t, n, d}{startRow:endRow, 5};  % total number of cells
                cells_fracDR{t, n, d, run} = dataCells{t, n, d}{startRow:endRow, 15} ./ dataCells{t, n, d}{startRow:endRow, 5};  % fraction of drug resistant cells
            end
        end
    end
end



% Calculate the mean for each heatmap point across runs
means_totcells = cell(3, length(doubling_time_resistant), length(dose_increments));
means_fracDR = cell(3, length(doubling_time_resistant), length(dose_increments));

for t = 1:3
    for n = 1:length(doubling_time_resistant_std)
        for d = 1:length(dose_increments)
            % Concatenate data for each cell and calculate the mean
            stack_totcells = cat(4, cells_totcells{t, n, d, :});
            stack_fracDR = cat(4, cells_fracDR{t, n, d, :});
            means_totcells{t, n, d} = mean(stack_totcells, 4);
            means_fracDR{t, n, d} = mean(stack_fracDR, 4);
        end
    end
end

%% 3. PLOT GRAPH

% Initialize matrices for heatmap data
data_for_heatmap = struct();

for t = 1:3
    data_for_heatmap.(dataTypes{t}) = struct(...
        'totcells', zeros(length(doubling_time_resistant), length(dose_increments)), ...
        'fracDR', zeros(length(doubling_time_resistant), length(dose_increments)) ...
        );
    
    for n = 1:length(doubling_time_resistant_std)
        for d = 1:length(dose_increments)
            data_for_heatmap.(dataTypes{t}).totcells(n, d) = means_totcells{t, n, d}(end);
            data_for_heatmap.(dataTypes{t}).fracDR(n, d) = means_fracDR{t, n, d}(end);
        end
    end
end

figure;

% Subplot 1: Single-cell Total Cells
subplot(2, 3, 1);
imagesc(rot90(data_for_heatmap.singlecell.totcells));
colormap('copper');
colorbar;
set(gca, 'XTick', 1:size(data_for_heatmap.singlecell.totcells, 2), ...
    'XTickLabel', {'41','49.2','57.4','65.6','73.8','82'}, ...
    'YTick', 1:size(data_for_heatmap.singlecell.totcells, 1), ...
    'YTickLabel', {'1','0.8','0.6','0.4','0.2','0'});
xlabel('Number of seeded clusters');
ylabel('Dose of drugs 1 and 2 (\muM)');
title('Single-cell Total Cell Count');

% Subplot 2: Multi-cell Total Cells
subplot(2, 3, 2);
imagesc(rot90(data_for_heatmap.multicell.totcells));
colormap('copper');
colorbar;
set(gca, 'XTick', 1:size(data_for_heatmap.multicell.totcells, 2), ...
    'XTickLabel', {'41','49.2','57.4','65.6','73.8','82'}, ...
    'YTick', 1:size(data_for_heatmap.multicell.totcells, 1), ...
    'YTickLabel', {'1','0.8','0.6','0.4','0.2','0'});
xlabel('Number of seeded clusters');
ylabel('Dose of drugs 1 and 2 (\muM)');
title('Multi-cell Total Cell Count');

% Subplot 3: Mono-clusters Total Cells
subplot(2, 3, 3);
imagesc(rot90(data_for_heatmap.monoclusters.totcells));
colormap('copper');
colorbar;
set(gca, 'XTick', 1:size(data_for_heatmap.monoclusters.totcells, 2), ...
    'XTickLabel', {'41','49.2','57.4','65.6','73.8','82'}, ...
    'YTick', 1:size(data_for_heatmap.monoclusters.totcells, 1), ...
    'YTickLabel', {'1','0.8','0.6','0.4','0.2','0'});
xlabel('Number of seeded clusters');
ylabel('Dose of drugs 1 and 2 (\muM)');
title('Mono-clusters Total Cell Count');

% Subplot 4: Single-cell Fraction DR
subplot(2, 3, 4);
imagesc(rot90(data_for_heatmap.singlecell.fracDR));
colormap('copper');
colorbar;
set(gca, 'XTick', 1:size(data_for_heatmap.singlecell.fracDR, 2), ...
    'XTickLabel', {'41','49.2','57.4','65.6','73.8','82'}, ...
    'YTick', 1:size(data_for_heatmap.singlecell.fracDR, 1), ...
    'YTickLabel', {'1','0.8','0.6','0.4','0.2','0'});
xlabel('Number of seeded clusters');
ylabel('Dose of drugs 1 and 2 (\muM)');
title('Single-cell Fraction of Drug Resistant Cells');

% Subplot 5: Multi-cell Fraction DR
subplot(2, 3, 5);
imagesc(rot90(data_for_heatmap.multicell.fracDR));
colormap('copper');
colorbar;
set(gca, 'XTick', 1:size(data_for_heatmap.multicell.fracDR, 2), ...
    'XTickLabel', {'41','49.2','57.4','65.6','73.8','82'}, ...
    'YTick', 1:size(data_for_heatmap.multicell.fracDR, 1), ...
    'YTickLabel', {'1','0.8','0.6','0.4','0.2','0'});
xlabel('Number of seeded clusters');
ylabel('Dose of drugs 1 and 2 (\muM)');
title('Multi-cell Fraction of Drug Resistant Cells');

% Subplot 6: Mono-clusters Fraction DR
subplot(2, 3, 6);
imagesc(rot90(data_for_heatmap.monoclusters.fracDR));
colormap('copper');
colorbar;
set(gca, 'XTick', 1:size(data_for_heatmap.monoclusters.fracDR, 2), ...
    'XTickLabel', {'41','49.2','57.4','65.6','73.8','82'}, ...
    'YTick', 1:size(data_for_heatmap.monoclusters.fracDR, 1), ...
    'YTickLabel', {'1','0.8','0.6','0.4','0.2','0'});
xlabel('Number of seeded clusters');
ylabel('Dose of drugs 1 and 2 (\muM)');
title('Mono-clusters Fraction of Drug Resistant Cells');