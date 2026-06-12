% =========================================================================
% SCRIPT: main_parcellation_pipeline_retest.m
% =========================================================================
% DESCRIPTION:
%   Batch parcellation pipeline for ALL HCP Retest Sessions.
%
%   Processes:
%       REST1_LR
%       REST1_RL
%       REST2_LR
%       REST2_RL
%
%   for all retest subjects simultaneously.
% =========================================================================

clc;
clear;
close all;

%% =========================================================================
% 1. GLOBAL CONFIGURATION
% =========================================================================

project_root = 'F:\PhD Code\My_PhD_Project';

wb_cmd_path = fullfile(project_root, ...
    'tools', 'workbench\bin_windows64\wb_command.exe');

atlas_name = 'Schaefer200_Kong17';

atlas_file = fullfile(project_root, ...
    'data', 'atlases', ...
    'Schaefer2018_200Parcels_Kong2022_17Networks_order.dlabel.nii');

% Retest raw data root
raw_data_root = 'D:\HCP\retest\Data';

% Add cifti toolbox
addpath(genpath(fullfile(project_root, ...
    'tools', 'cifti-matlab-master')));

%% =========================================================================
% 2. DEFINE ALL SESSIONS
% =========================================================================

sessions = { ...
    'REST1_LR', ...
    'REST1_RL', ...
    'REST2_LR', ...
    'REST2_RL'};

%% =========================================================================
% 3. GET SUBJECT LIST
% =========================================================================

dir_content = dir(raw_data_root);

subjects = dir_content( ...
    [dir_content.isdir] & ...
    ~ismember({dir_content.name}, {'.', '..'}));

fprintf('=====================================================\n');
fprintf('HCP RETEST PARCELLATION PIPELINE\n');
fprintf('Subjects Found : %d\n', length(subjects));
fprintf('Sessions       : %d\n', length(sessions));
fprintf('=====================================================\n');

t_total = tic;

%% =========================================================================
% 4. MAIN LOOP
% =========================================================================

for s = 1:length(sessions)

    CURRENT_SESSION = sessions{s};

    fprintf('\n=====================================================\n');
    fprintf('SESSION: %s\n', CURRENT_SESSION);
    fprintf('=====================================================\n');

    % -------------------------------------------------
    % Dynamic Session Definitions
    % -------------------------------------------------

    internal_dir = ['rfMRI_' CURRENT_SESSION];

    file_suffix = sprintf( ...
        'rfMRI_%s_Atlas_MSMAll_hp2000_clean.dtseries.nii', ...
        CURRENT_SESSION);

    % -------------------------------------------------
    % Output Directory
    % -------------------------------------------------

    output_dir = fullfile( ...
        project_root, ...
        'data', ...
        'inter_parcellated', ...
        atlas_name, ...
        'retest', ...
        CURRENT_SESSION);

    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % -------------------------------------------------
    % Subject Loop
    % -------------------------------------------------

    count_ok = 0;

    for i = 1:length(subjects)

        subID = subjects(i).name;

        % Input dtseries
        input_nii = fullfile( ...
            raw_data_root, ...
            subID, ...
            internal_dir, ...
            file_suffix);

        % Output MAT
        output_mat = fullfile( ...
            output_dir, ...
            [subID '.mat']);

        % -------------------------------------------------
        % Check input existence
        % -------------------------------------------------

        if ~exist(input_nii, 'file')

            fprintf('[MISSING] %s | %s\n', ...
                subID, CURRENT_SESSION);

            continue;
        end

        % -------------------------------------------------
        % Resume capability
        % -------------------------------------------------

        if exist(output_mat, 'file')

            fprintf('[SKIP] %s | %s\n', ...
                subID, CURRENT_SESSION);

            continue;
        end

        % -------------------------------------------------
        % Run Parcellation
        % -------------------------------------------------

        fprintf('Processing %s | %s ... ', ...
            subID, CURRENT_SESSION);

        [is_success, msg] = parcellate_cifti_wrapper( ...
            input_nii, ...
            atlas_file, ...
            output_mat, ...
            wb_cmd_path);

        if is_success

            fprintf('DONE\n');
            count_ok = count_ok + 1;

        else

            fprintf('FAILED (%s)\n', msg);

        end

    end

    fprintf('\nSession %s Finished | Success: %d\n', ...
        CURRENT_SESSION, count_ok);

end

fprintf('\n=====================================================\n');
fprintf('ALL RETEST SESSIONS FINISHED\n');
fprintf('Total Time: %.2f mins\n', toc(t_total)/60);
fprintf('=====================================================\n');