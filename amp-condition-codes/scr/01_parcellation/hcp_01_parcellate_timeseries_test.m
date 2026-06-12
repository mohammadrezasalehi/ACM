% =========================================================================
% SCRIPT: hcp_01_parcellate_timeseries.m
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

%% 1. CONFIGURATION
% =========================================================================
% Define Base Paths (Use absolute paths or relative to script location)
project_root = 'F:\PhD Code\My_PhD_Project'; % Adjust this to your root

% External Tool Paths
path_wb_command = fullfile(project_root, 'tools', 'workbench\bin_windows64\wb_command.exe');
path_cifti_toolbox = fullfile(project_root, 'tools', 'cifti-matlab-master');

% Input Data Settings
data_dir_name = 'R1LR'; % Example: 'R1LR', 'R1RL'
input_data_root = 'D:\HCP\retest\Data'; % Path to raw HCP subject folders

% Atlas Settings
atlas_file = fullfile(project_root, 'data', 'atlases', ...
    'Schaefer2018_200Parcels_Kong2022_17Networks_order.dlabel.nii');

% Output Settings
output_dir = fullfile(project_root, 'data', 'inter_parcellated', 'kong200', data_dir_name);
temp_dir = fullfile(output_dir, 'temp_cifti'); % For intermediate .ptseries files

% Setup Environment
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
if ~exist(temp_dir, 'dir'), mkdir(temp_dir); end
addpath(genpath(path_cifti_toolbox));

%% 2. SUBJECT DISCOVERY
% =========================================================================
fprintf('Scanning for subjects in: %s\n', input_data_root);
dir_content = dir(input_data_root);
% Filter only directories and exclude '.' and '..'
subject_list = dir_content([dir_content.isdir] & ~ismember({dir_content.name}, {'.', '..'}));

fprintf('Found %d subjects.\n', length(subject_list));

%% 3. MAIN PROCESSING LOOP
% =========================================================================
tic;
counter = 0;

for i = 1:length(subject_list)
    subID = subject_list(i).name;
    
    % Construct file paths
    % Note: Adjust the internal folder structure ('rfMRI_REST1_LR' etc.) as needed
    nii_filename = 'rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii';
    input_path = fullfile(input_data_root, subID, 'rfMRI_REST1_LR', nii_filename);
    
    output_ptseries_path = fullfile(temp_dir, [subID, '_parcellated.ptseries.nii']);
    output_mat_path = fullfile(output_dir, [subID, '.mat']);
    
    % Check if input exists and output doesn't (to allow resume)
    if exist(input_path, 'file') && ~exist(output_mat_path, 'file')
        
        fprintf('Processing Subject: %s ... ', subID);
        
        try
            % -------------------------------------------------------------
            % A. Run wb_command -cifti-parcellate
            % -------------------------------------------------------------
            % Using "COLUMN" direction for standard dense timeseries
            cmd = sprintf('"%s" -cifti-parcellate "%s" "%s" COLUMN "%s"', ...
                path_wb_command, input_path, atlas_file, output_ptseries_path);
            
            [status, cmdout] = system(cmd);
            
            if status ~= 0
                error('wb_command failed: %s', cmdout);
            end
            
            % -------------------------------------------------------------
            % B. Load Parcellated Data & Save as .mat
            % -------------------------------------------------------------
            cifti_data = ft_read_cifti(output_ptseries_path);
            
            % Extract Data (Rows: Parcels, Cols: Time)
            % Depending on the atlas and cifti-matlab version, data might be 
            % in .ptseries or .cdata.
            if isfield(cifti_data, 'ptseries')
                data_on_atlas = cifti_data.ptseries;
            elseif isfield(cifti_data, 'cdata')
                data_on_atlas = cifti_data.cdata;
            else
                error('Could not find data field in loaded cifti structure.');
            end
            
            % Save
            save(output_mat_path, 'data_on_atlas');
            
            % Cleanup temp file
            delete(output_ptseries_path);
            
            fprintf('Done.\n');
            counter = counter + 1;
            
        catch ME
            fprintf('\n[ERROR] Failed on %s: %s\n', subID, ME.message);
        end
        
    elseif exist(output_mat_path, 'file')
        fprintf('Skipping %s (Already processed).\n', subID);
    else
        % Input file missing silently (or log it if needed)
        % fprintf('Input missing for %s\n', subID);
    end
end

%% 4. SUMMARY
% =========================================================================
elapsed_time = toc;
fprintf('\nProcessing Completed.\n');
fprintf('Total Subjects Processed: %d\n', counter);
fprintf('Total Time: %.2f minutes\n', elapsed_time/60);