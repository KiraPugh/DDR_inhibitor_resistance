load('output_abm_calibration');

% List of variable names
varNames = {'dmso', 'cera_01', 'cera_03', 'cera_1', 'ola_02', 'ola_1', ...
    'comb_01_1', 'comb_03_02', 'comb_03_1', 'comb_1_02', 'comb_1_1'};

time = (0:1:310)';

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

for i = 1:6
    subplot(2, 3, i);
    
    for j = 1:size(data_sets{i}, 1)
        varName = data_sets{i}{j, 1};
        color = data_sets{i}{j, 2};
        
        % Create dynamic variable names for mean and std
        mean_var_name = sprintf('mean_%s', varName);
        std_var_name = sprintf('std_%s', varName);
        
        % Retrieve the data using eval to access the dynamically named variables
        mean_data = eval(mean_var_name);
        std_data = eval(std_var_name);
        
        mean_numeric = cell2mat(mean_data);
        std_numeric = cell2mat(std_data);
        
        fill([time; flipud(time)], [mean_numeric + std_numeric; flipud(mean_numeric - std_numeric)], ...
            color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        hold on;
        
        % Plot the mean line
        plot(time, mean_numeric, 'Color', color, 'LineWidth', 2);
        hold on;
        
    end
    
    xlabel('Time (hours)');
    ylabel('Cell Confluency');
    set(gca, 'FontSize', 15);
    title(titles{i});
    xlim([0, 310]);
end
