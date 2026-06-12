% Generates Surrogate Multivariate Gaussian Data preserving Static FC.
% Runs the dynamic modulation metrics on the surrogate data.

clc; clear; close all;

%% 1. Configuration
further_cleaning = 'inter_parcellated'; % Adjust based on your pipeline
session = 'REST1_LR';
cleaning_GSR = 'With_GSR\'; % or 'No_clean\'

% Input Data Directory
data_root = fullfile('F:\PhD Code\My_PhD_Project\data', further_cleaning, ...
        'Schaefer200_Kong17\', session); % Adjust if needed

% NEW Output Directory for Gaussian Null
output_dir = fullfile('F:\PhD Code\My_PhD_Project\results\Modulated_FC_Gaussian\', ...
        session, 'No_clean\');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

subjects = dir(fullfile(data_root, '*.mat')); 
fprintf('Starting Gaussian Null Generation for %d subjects...\n', length(subjects));

%% 2. Processing Loop
h = waitbar(0, 'Generating Gaussian Surrogates & Running Metrics...');

for i = 1:length(subjects)
    subj_name = subjects(i).name;
    subj_path = fullfile(subjects(i).folder, subj_name);
    
    % A. Load Real Data
    tmp = load(subj_path);
    ts_real = tmp.data_on_atlas'; % [1200 Timepoints x 200 Nodes]
    [ntime, nnodes] = size(ts_real);
    
    % B. Generate Multivariate Gaussian Surrogate
    % 1. Calculate Empirical Mean and Covariance
    mu = mean(ts_real, 1);
    Sigma = cov(ts_real);
    
    % Ensure Sigma is perfectly symmetric to avoid mvnrnd floating point errors
    Sigma = (Sigma + Sigma') / 2; 
    
    % 2. Generate random data from this distribution
    % rng('shuffle'); % Optional: ensure different random seeds
    ts_gaussian = mvnrnd(mu, Sigma, ntime);
    
    % C. Run Your Metric on Gaussian Data
    % (Using the optimized ultra-fast version we wrote earlier)
    subj_metrics = calc_dynamic_modulation_metrics(ts_gaussian, 'n_levels', 20); 
    
    % D. Save Results
    [~, name_no_ext, ~] = fileparts(subj_name);
    save_path = fullfile(output_dir, [name_no_ext '_dyn_mod.mat']);
    save(save_path, '-struct', 'subj_metrics');
    
    waitbar(i / length(subjects), h, sprintf('Processed %d/%d', i, length(subjects)));
end
close(h);
fprintf('Gaussian Null Batch Complete!\n');