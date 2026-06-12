% =========================================================================
% SCRIPT: main_parcellation_pipeline.m
% =========================================================================
% DESCRIPTION:
%   Master script to batch process HCP data across different sessions.
%   Controls the 'parcellate_cifti_wrapper' function.
% =========================================================================
% AUTHOR:       Mohammad Reza Salehi
% DATE:         2025-11-26]
% VERSION:      1.0
% MATLAB VER:   R2018b
%
% DESCRIPTION:
%   This script iterates through HCP rfMRI .dtseries.nii files and 
%   parcellates them onto a specified atlas using Connectome Workbench 
%   (wb_command). The output is then converted to a .mat file containing
%   the average time series for each parcel.
%
% INPUTS:
%   - HCP "clean" dtseries files (Surface-mapped MSMAll)
%   - Schaefer 2018 (Kong 2022 Networks) Atlas (.dlabel.nii)
%
% OUTPUTS:
%   - .mat files containing (N_Parcels x TimePoints) matrix for each subject.
%
% DEPENDENCIES:
%   1. Connectome Workbench (wb_command) installed.
%   2. cifti-matlab toolbox (or FieldTrip) in MATLAB path.
%
% USAGE:
%   Update the "CONFIGURATION" section below with your local paths before
%   running.
% =========================================================================

clc;
clear;
close all;

%% 1. GLOBAL CONFIGURATION (USER EDIT HERE)
% =========================================================================

% ---> SELECT THE SESSION TO PROCESS <---
% Options: 'REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL'
CURRENT_SESSION = 'REST1_LR'; 

% Project Roots
project_root = 'F:\PhD Code\My_PhD_Project'; % Adjust to your path
wb_cmd_path = fullfile(project_root, 'tools', 'workbench\bin_windows64\wb_command.exe');

% Atlas Settings
atlas_name = 'Schaefer200_Kong17'; % Used for output folder naming
atlas_file = fullfile(project_root, 'data', 'atlases', ...
    'Schaefer2018_200Parcels_Kong2022_17Networks_order.dlabel.nii');

% Toolbox
addpath(genpath(fullfile(project_root, 'tools', 'cifti-matlab-master')));


%% 2. SESSION MAPPING (Dynamic Path Generation)
% =========================================================================
% Here we define where the data lives for each session.
% You can add more sessions here easily.

switch CURRENT_SESSION
    case 'REST1_LR'
        raw_data_root = 'D:\HCP\retest\Data';     % Where folders like 100307 live
        internal_dir  = 'rfMRI_REST1_LR';         % Folder inside subject dir
        file_suffix   = 'rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii';
        
    case 'REST1_RL'
        raw_data_root = 'D:\HCP\R1RL\Data'; 
        internal_dir  = 'rfMRI_REST1_RL';
        file_suffix   = 'rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii';
        
    case 'REST2_LR'
        raw_data_root = 'D:\HCP\R2LR\Data';
        internal_dir  = 'rfMRI_REST2_LR';
        file_suffix   = 'rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii';
        
    case 'REST2_RL'
        raw_data_root = 'D:\HCP\R2RL\Data';
        internal_dir  = 'rfMRI_REST2_RL';
        file_suffix   = 'rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii';
        
    otherwise
        error('Session %s not defined in switch case.', CURRENT_SESSION);
end

% Construct Output Directory Automatically
% Structure: data/inter_parcellated/<AtlasName>/<SessionName>/
output_dir = fullfile(project_root, 'data', 'inter_parcellated', atlas_name, CURRENT_SESSION);

if ~exist(output_dir, 'dir'), mkdir(output_dir); end


%% 3. BATCH PROCESSING LOOP
% =========================================================================
fprintf('======================================================\n');
fprintf('STARTING PIPELINE\n');
fprintf('Session: %s\n', CURRENT_SESSION);
fprintf('Input:   %s\n', raw_data_root);
fprintf('Output:  %s\n', output_dir);
fprintf('======================================================\n');

% Get subject list
dir_content = dir(raw_data_root);
subjects = dir_content([dir_content.isdir] & ~ismember({dir_content.name}, {'.', '..'}));

t_start = tic;
count_ok = 0;

for i = 1:length(subjects)
    subID = subjects(i).name;
    
    % Define File Paths
    input_nii = fullfile(raw_data_root, subID, internal_dir, file_suffix);
    output_mat = fullfile(output_dir, [subID, '.mat']);
    
    % Resume Capability: Skip if exists
    if exist(output_mat, 'file')
        % fprintf('Skipping %s (Exists)\n', subID);
        continue;
    end
    
    fprintf('Processing %s ... ', subID);
    
    % CALL THE WRAPPER FUNCTION
    [is_success, msg] = parcellate_cifti_wrapper(input_nii, atlas_file, output_mat, wb_cmd_path);
    
    if is_success
        fprintf('DONE.\n');
        count_ok = count_ok + 1;
    else
        fprintf('FAILED. (%s)\n', msg);
    end
end

fprintf('\nPipeline Finished. Processed %d subjects in %.2f mins.\n', count_ok, toc(t_start)/60);