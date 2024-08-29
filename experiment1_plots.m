%% 1. IMPORT THE DATA

% Define the parameters
DRfrac_values = [10, 30];
DRtype = 3;
DoseC_values = [0, 30, 100];
DoseO_values = [0, 20, 100];
numClusters_values = [419, 32, 1]; % For single-cell clusters, multi-cell clusters, and monoclusters
numRuns = 100;
rowsPerRun = 311;

% Create the cell arrays to hold the data
combres_cells = struct();

% Loop through each combination of parameters
for DRfrac = DRfrac_values
    for numClusters = numClusters_values
        for DoseC = DoseC_values
            for DoseO = DoseO_values
                
                % Skip the combination where DoseC = 30, DoseO = 100, and numClusters = 419
                
                if (numClusters == 419 || numClusters == 32 || numClusters == 1) && ...
                        (DoseC == 30 && DoseO == 100 || DoseC == 100 && DoseO == 20)
                    continue; % Skip this iteration
                end
                
                
                % Construct the file name
                fileName = sprintf('oct_B_mu41000_sigma8200_DoseC%d_DoseO%d_NoCircles%d_DRfrac%d_DRtype%d', ...
                    DoseC, DoseO, numClusters, DRfrac, DRtype);
                %%%combres_10s_dmso_cells = cell(numRuns, 1);
                
                
                % Generate the field name string based on the parameters
                fieldName = sprintf('combres_%d_%d_%dc_%do', DRfrac, numClusters, DoseC, DoseO);
                
                % Load the table from the file
                dataTable = readtable(fileName);
                
                % Store the table as a variable in the base workspace
                assignin('base', fieldName, dataTable);
                
                fieldName2 = sprintf('combres_cells_%d_%d_%dc_%do', DRfrac, numClusters, DoseC, DoseO);
                fieldName3 = sprintf('combres_cells2_%d_%d_%dc_%do', DRfrac, numClusters, DoseC, DoseO);
                
                fieldName_stack = sprintf('combres_stack_%d_%d_%dc_%do', DRfrac, numClusters, DoseC, DoseO);
                fieldName_mean = sprintf('combres_mean_%d_%d_%dc_%do', DRfrac, numClusters, DoseC, DoseO);
                fieldName_std = sprintf('combres_std_%d_%d_%dc_%do', DRfrac, numClusters, DoseC, DoseO);
                
                % Create the cell array and assign it to the base workspace
                % Create the cell array and assign it to the base workspace
                dataCell = cell(numRuns, 1);
                dataCell2 = cell(numRuns, 1);
                
                % Split the dataTable into `numRuns` cells and perform calculations
                for run = 1:numRuns
                    startRow = (run - 1) * rowsPerRun + 1;
                    endRow = run * rowsPerRun;
                    
                    % Perform the calculation directly when filling the cell array
                    % Ensure that columns 15 and 5 exist in the dataTable
                    if size(dataTable, 2) >= 15
                        % Perform the calculation and store the result
                        dataCell{run} = dataTable{startRow:endRow, 15} ./ dataTable{startRow:endRow, 5} * 100;
                        dataCell2{run} = dataTable{startRow:endRow, 13} ./ dataTable{startRow:endRow, 5} * 100;
                        
                        dataCell_stack = cat(3,dataCell{:}); % Stack data from all runs along the third dimension
                        dataCell_mean = mean(dataCell_stack, 3);
                        dataCell_std = std(dataCell_stack, 0, 3);
                        
                        dataCell_stack2 = cat(3,dataCell2{:}); % Stack data from all runs along the third dimension
                        dataCell_mean2 = mean(dataCell_stack2, 3);
                        dataCell_std2 = std(dataCell_stack2, 0, 3);
                    else
                        % Handle the case where there are fewer than 15 columns
                        error('dataTable does not have enough columns.');
                    end
                end
                
                % Assign the cell array to the base workspace
                assignin('base', fieldName2, dataCell);
                assignin('base', fieldName3, dataCell2);
                
                assignin('base', fieldName_stack, dataCell_stack);
                assignin('base', fieldName_mean, dataCell_mean);
                assignin('base', fieldName_std, dataCell_std);
                
                
                % Plotting
                figure;
                
                % Ensure 'time' is properly defined
                time = [0:1:310]';
                
                % Plot the mean and standard deviation
                fill_area = fill([time; flipud(time)], [dataCell_mean + dataCell_std; flipud(dataCell_mean - dataCell_std)], 'w');
                set(fill_area, 'EdgeColor', 'none'); % Remove the edge line
                alpha(fill_area, 0.3); % Adjust transparency
                
                hold on;
                
                h = area(time, [dataCell_mean, dataCell_mean2], 'LineWidth', 1.5);
    
    set(h(1), 'FaceColor', [153/255, 142/255, 195/255]); 
    set(h(2), 'FaceColor', [241/255, 163/255, 64/255]);
    
    fill_area = fill([time; flipud(time)], [dataCell_mean + dataCell_std; flipud(dataCell_mean + dataCell_std)], 'w');
    set(fill_area, 'EdgeColor', 'none'); % Remove the line on the fill area
    alpha(fill_area, 0.3); % Adjust the alpha value here
    
                %plot(time, dataCell_mean, 'LineWidth', 1.5, 'Color', [0 0 1]); % Mean plot
                plot(time, dataCell_mean + dataCell_std, 'k--', 'LineWidth', 1); % Mean + std plot
                plot(time, dataCell_mean - dataCell_std, 'k--', 'LineWidth', 1); % Mean - std plot
                
                
                % Add labels, title, etc.
                xlabel('Time (hours)');
                ylabel('Fraction');
                title(sprintf('DRfrac %d, NumClusters %d, DoseC %d, DoseO %d', DRfrac, numClusters, DoseC, DoseO));
                grid on;
                xlim([0 310]);
                ylim([0 100]);
                
                ax = gca;
                ax.YAxis.Exponent = 0;
                set(gca,'XTick',[],'YTick',[]);
                
                hold off;
                
                
                
                
            end
        end
    end
end





