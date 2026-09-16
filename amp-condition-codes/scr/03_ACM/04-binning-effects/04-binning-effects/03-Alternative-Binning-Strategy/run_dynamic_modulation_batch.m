%% Suggested Name: run_dynamic_modulation_batch.m
% Batch script to calculate Amplitude-based Dynamic Modulation Metrics for all subjects.

clc; clear; close all;

%% 1. Configuration
% Select Session and Preprocessing Pipline
% Valid for Sessions: REST1_LR, REST1_RL, REST2_LR, REST2_RL
% Valid for Preprocessing: No_GSR, With_GSR
% Valid for further_cleaning: processed_clean, inter_parcellated

further_cleaning = 'inter_parcellated';
session = 'REST1_LR';
cleaning_GSR = 'With_GSR\';

% --- NEW SETTINGS FOR ALTERNATIVE BINNING STRATEGIES ---
% Choose the strategy: 'Quantile' (Original), 'Uniform', 'SlidingWindow'
bin_strategy = 'SlidingWindow'; % Change to 'SlidingWindow'
bin_number = 20; % Keep 20 to compare fairly with the original 20

% Dynamic Folder and File Naming
folder_suffix = sprintf('Strategy_%s', bin_strategy); 
file_suffix = sprintf('_dyn_mod_%s.mat', bin_strategy); 

% Set the base directory dynamically
base_result_dir = fullfile('F:\PhD Code\My_PhD_Project\results\Modulated_FC_Alternatives', folder_suffix);

if strcmp(further_cleaning, 'processed_clean') 
    data_root = fullfile('F:\PhD Code\My_PhD_Project\data', further_cleaning, ...
        'Schaefer200_Kong17\', session, cleaning_GSR);
    output_dir = fullfile(base_result_dir, session, cleaning_GSR);
    
elseif strcmp(further_cleaning, 'inter_parcellated')
    data_root = fullfile('F:\PhD Code\My_PhD_Project\data', further_cleaning, ...
        'Schaefer200_Kong17\', session);
    output_dir = fullfile(base_result_dir, session, 'No_clean\');
end

% Get list of files
subjects = dir(fullfile(data_root, '*.mat')); 

% Create output directory if it doesn't exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir); 
    fprintf('Created new directory: %s\n', output_dir);
end

%% 2. Analysis Loop
fprintf('Starting Dynamic Modulation batch analysis for %d subjects...\n', length(subjects));

% Start Waitbar
h = waitbar(0, 'Initializing Dynamic Modulation Analysis...');

for i = 1:length(subjects)
    subj_name = subjects(i).name;
    subj_path = fullfile(subjects(i).folder, subj_name);
    
    % Load Data
    tmp = load(subj_path);
    
    if strcmp(further_cleaning, 'processed_clean') 
        ts = tmp.clean_data'; % Adjust based on your file structure
    elseif strcmp(further_cleaning, 'inter_parcellated')
        ts = tmp.data_on_atlas'; % Adjust based on your file structure
    end
    
    
    % --- CORE CALCULATION ---
    if strcmp(bin_strategy, 'Quantile')
        subj_metrics = calc_dynamic_modulation_metrics(ts, 'n_levels', bin_number); 
    elseif strcmp(bin_strategy, 'Uniform')
        subj_metrics = calc_uniform_amplitude_metrics(ts, 'n_levels', bin_number);
        fprintf('   -> Average Empty Bins (NaNs) created by Uniform Strategy: %.2f%%\n',...
            mean(subj_metrics.nan_percent_per_node));
    elseif strcmp(bin_strategy, 'SlidingWindow')
        % Call the new function. 
        % Window = 120 (same stability as K=10), Step = 20 (Yields 55 overlapping points)
        subj_metrics = calc_sliding_window_metrics(ts, 'window_size', 60, 'step_size', 10); 
    end
    
    % Save Results
    [~, name_no_ext, ~] = fileparts(subj_name);
    save_path = fullfile(output_dir, [name_no_ext file_suffix]);
    save(save_path, '-struct', 'subj_metrics');
    
    % Update Waitbar
    waitbar(i / length(subjects), h, sprintf('Processed %d/%d: %s', i, length(subjects), name_no_ext));
end

close(h);
fprintf('Batch analysis complete. Results successfully saved in:\n%s\n', output_dir);