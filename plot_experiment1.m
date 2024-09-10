%% IMPORT THE DATA

%load('output_exp1_ATRiPARPires.mat'); % Plot Fig 3
%load('output_exp1_ATRires.mat'); % Plot Fig 4
load('output_exp1_PARPires.mat'); % Plot Fig 5


time = [0:1:310]';

% Define the conditions and corresponding data arrays
conditions = {'10s_dmso', '10mc_dmso', '10c_dmso', '30s_dmso', '30mc_dmso', '30c_dmso', ...
    '10s_cera_01', '10s_cera_03', '10s_ola_003', '10s_ola_03', '10s_comb_01_003', '10s_comb_03_03', ...
    '30s_cera_01', '30s_cera_03', '30s_ola_003', '30s_ola_03', '30s_comb_01_003', '30s_comb_03_03', ...
    '10mc_cera_01', '10mc_cera_03', '10mc_ola_003', '10mc_ola_03', '10mc_comb_01_003', '10mc_comb_03_03', ...
    '30mc_cera_01', '30mc_cera_03', '30mc_ola_003', '30mc_ola_03', '30mc_comb_01_003', '30mc_comb_03_03', ...
    '10c_cera_01', '10c_cera_03', '10c_ola_003', '10c_ola_03', '10c_comb_01_003', '10c_comb_03_03', ...
    '30c_cera_01', '30c_cera_03', '30c_ola_003', '30c_ola_03', '30c_comb_01_003', '30c_comb_03_03'};

num_conditions = length(conditions);

%% PLOT THE FRACTION OF DRUG RESISTANT CELLS

% Number of rows and columns for the subplot grid
rows = 7;
cols = 6;

figure;

for i = 1:num_conditions
    
    subplot(rows, cols, i);
    
    % Create the variable names dynamically
    var_name_1 = ['mean_combres_' conditions{i}];
    var_name_2 = ['mean_combres_' conditions{i} '_2'];
    
    % Evaluate the expressions
    data_1 = eval(var_name_1);
    data_2 = eval(var_name_2);
    
    % Create a filled region between data_1 and data_2
    h = area(time, [data_1, data_2], 'LineWidth', 1.5);
    
    % Set colors for the areas
    set(h(1), 'FaceColor', [153/255, 142/255, 195/255]);
    set(h(2), 'FaceColor', [241/255, 163/255, 64/255]);
    
    hold on;
    
    % Plot the shaded error region around data_1
    std_data = eval(['std_combres_' conditions{i}]);
    fill_area = fill([time; flipud(time)], ...
        [data_1 + std_data; flipud(data_1 - std_data)], ...
        'w', 'EdgeColor', 'none');
    alpha(fill_area, 0.3);
    
    % Overlay dashed lines for the error margins
    plot(time, data_1 + std_data, 'k--', 'LineWidth', 1);
    plot(time, data_1 - std_data, 'k--', 'LineWidth', 1);
    
    % Remove grid and tick marks
    set(gca, 'XTick', [], 'YTick', []);
    
    % Set axis limits
    xlim([0 310]);
    ylim([0 100]);
    
    title([conditions{i}], 'Interpreter', 'none');
    
    hold off;
end


%% PLOT THE TOTAL NUMBER OF CELLS

% Number of rows and columns in the subplot grid
rows = 7;
cols = 6;

figure;

for i = 1:num_conditions
    
    subplot(rows, cols, i);
    
    % Construct the variable names
    std_var_name = ['std_combres_' conditions{i} '_3'];
    mean_var_name = ['mean_combres_' conditions{i} '_3'];
    
    % Evaluate the variables
    std_totalcells = eval(std_var_name);
    mean_totalcells = eval(mean_var_name);
    
    % Plot the shaded error region
    fill([time; flipud(time)], [mean_totalcells + std_totalcells; flipud(mean_totalcells - std_totalcells)], ...
        'r', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    hold on;
    
    % Plot the mean line
    plot(time, mean_totalcells, 'Color', 'r', 'LineWidth', 2);
    hold on;
    
    xlim([0 310]);
    
    title(conditions{i}, 'Interpreter', 'none');
end
