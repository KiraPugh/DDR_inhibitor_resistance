
%% 1. IMPORT THE DATA

combres_10s_dmso = readtable('oct_B_mu41000_sigma8200_DoseC0_DoseO0_NoCircles419_DRfrac30_DRtype3');
combres_10s_cera_01 = readtable('oct_B_mu41000_sigma8200_DoseC10_DoseO0_NoCircles419_DRfrac0_DRtype0');
combres_10s_cera_03 = readtable('oct_B_mu41000_sigma8200_DoseC30_DoseO0_NoCircles419_DRfrac0_DRtype0');
combres_10s_cera_1 = readtable('oct_B_mu41000_sigma8200_DoseC100_DoseO0_NoCircles419_DRfrac0_DRtype0');
combres_10s_ola_02 = readtable('oct_B_mu41000_sigma8200_DoseC0_DoseO20_NoCircles419_DRfrac0_DRtype0');
combres_10s_ola_1 = readtable('oct_B_mu41000_sigma8200_DoseC0_DoseO100_NoCircles419_DRfrac0_DRtype0');
combres_10s_comb_01_1 = readtable('oct_B_mu41000_sigma8200_DoseC10_DoseO100_NoCircles419_DRfrac0_DRtype0');
combres_10s_comb_03_02 = readtable('oct_B_mu41000_sigma8200_DoseC30_DoseO20_NoCircles419_DRfrac0_DRtype0');
combres_10s_comb_03_1 = readtable('oct_B_mu41000_sigma8200_DoseC30_DoseO100_NoCircles419_DRfrac0_DRtype0');
combres_10s_comb_1_02 = readtable('oct_B_mu41000_sigma8200_DoseC100_DoseO20_NoCircles419_DRfrac0_DRtype0');
combres_10s_comb_1_1 = readtable('oct_B_mu41000_sigma8200_DoseC100_DoseO100_NoCircles419_DRfrac0_DRtype0');




%% 2. SORT OUT THE INDIVIDUAL RUNS OF THE DATA

% Define the number of runs and the number of rows per run
numRuns = 3;
rowsPerRun = 311;


%% Comb Res 
% List of variable names
% List of variable names
varNames = {'dmso', 'cera_01', 'cera_03', 'cera_1', 'ola_02', 'ola_1', ...
            'comb_01_1', 'comb_03_02', 'comb_03_1', 'comb_1_02', 'comb_1_1'};

% Initialize cell arrays
for i = 1:length(varNames)
    eval(sprintf('combres_10s_%s_cells = cell(numRuns, 1);', varNames{i}));
end

% Loop through the runs
for run = 1:numRuns
    startRow = (run - 1) * rowsPerRun + 1;
    endRow = run * rowsPerRun;
    
    for i = 1:length(varNames)
        varName = varNames{i};
        eval(sprintf('combres_10s_%s_cells{run} = combres_10s_%s{startRow:endRow, 6}/100;', varName, varName));
    end
end

% Initialize struct to store means and standard deviations
mean_combres_10s = struct();
std_combres_10s = struct();

% Calculate means and standard deviations
for i = 1:length(varNames)
    varName = varNames{i};
    % Stack data from all runs along the third dimension
    eval(sprintf('combres_10s_%s_stack = cat(3, combres_10s_%s_cells{:});', varName, varName));
    
    % Calculate the mean and standard deviation
    eval(sprintf('mean_combres_10s.%s = mean(combres_10s_%s_stack, 3);', varName, varName));
    eval(sprintf('std_combres_10s.%s = std(combres_10s_%s_stack, 0, 3);', varName, varName));
end






% Time vector
time = [0:1:310]';

% Define subplot titles and data sets
titles = {'(I) Cell Confluency - DMSO', ...
          '(II) Cell Confluency - PARPi', ...
          '(III) Cell Confluency - ATRi', ...
          '(IV) Cell Confluency - Combinations', ...
          '(V) Cell Confluency - Combinations', ...
          '(VI) Cell Confluency - Combinations'};

data_sets = { {'dmso', 'k'}, ...
              {'ola_02', 'r'; ...
               'ola_1', [1, 0.5, 0]}, ...
              {'cera_01', [0.49, 0.18, 0.56]; ...
               'cera_03', 'g'; ...
               'cera_1', 'c'}, ...
              {'comb_01_1', 'b'}, ...
              {'comb_03_02', [0.93, 0.69, 0.13]; ...
               'comb_03_1', [0.72, 0.27, 1.00]}, ...
              {'comb_1_1', [0.64, 0.08, 0.18]; ...
               'comb_1_02', [0.47, 0.67, 0.19]}};

figure;

% Loop through subplots
for i = 1:6
    subplot(2, 3, i);
    
    % Loop through each dataset in the current subplot
    for j = 1:size(data_sets{i}, 1)
        varName = data_sets{i}{j, 1};
        color = data_sets{i}{j, 2};
        
        % Get the mean and std from the struct
        mean_data = mean_combres_10s.(varName);
        std_data = std_combres_10s.(varName);
        
        % Plot the shaded error region
        fill([time; flipud(time)], [mean_data + std_data; flipud(mean_data - std_data)], ...
             color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        hold on;
        
        % Plot the mean line
        plot(time, mean_data, 'Color', color, 'LineWidth', 2);
        hold on;
        
    end
    
    xlabel('Time (hours)');
    ylabel('Cell Confluency');
    set(gca, 'FontSize', 15);
    title(titles{i});
    xlim([0, 310]);
end
