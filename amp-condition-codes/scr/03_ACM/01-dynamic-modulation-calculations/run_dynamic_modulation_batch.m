% Batch script to calculate Amplitude-based Dynamic Modulation Metrics for all subjects.

clc; clear; close all;

%% 1. Configuration
% Select Session and Preprocessing Pipline
% Valid for Sessions: REST1_LR, REST1_RL, REST2_LR, REST2_RL
% Valid for Preprocessing: No_GSR, With_GSR
% Valid for further_cleaning: processed_clean, inter_parcellated

further_cleaning = 'inter_parcellated';
session = 'REST1_RL';
cleaning_GSR = 'With_GSR\';

% Base directory for results of this specific analysis
base_result_dir = 'F:\PhD Code\My_PhD_Project\results\Modulated_FC\';

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
    % Calculate dynamic metrics (variance, slope, zcr, best_poly)
    % The function uses 'parfor' internally for speed.
    subj_metrics = calc_dynamic_modulation_metrics(ts, 'n_levels', 20, ...
        'req_metrics', {'all'}); 
    
    % --- Save Results ---
    [~, name_no_ext, ~] = fileparts(subj_name);
    % Added '_dyn_mod' suffix to easily identify the output files
    save_path = fullfile(output_dir, [name_no_ext '_dyn_mod.mat']);
    
    % Saving struct directly so variables are loaded neatly later
    save(save_path, '-struct', 'subj_metrics');
    
    % Update Waitbar
    waitbar(i / length(subjects), h, sprintf('Processed %d/%d: %s', i, length(subjects), name_no_ext));
end

close(h);
fprintf('Batch analysis complete. Results successfully saved in:\n%s\n', output_dir);