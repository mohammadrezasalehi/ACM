%% =========================================================================
% SCRIPT: run_dynamic_modulation_batch_retest.m
% =========================================================================
% DESCRIPTION:
%   Batch extraction of dynamic modulation metrics for ALL HCP RETEST data.
%
%   Processes all sessions together:
%       REST1_LR
%       REST1_RL
%       REST2_LR
%       REST2_RL
%
% OUTPUT STRUCTURE:
%
% results/Modulated_FC/retest/<SESSION>/<PIPELINE>/<SUBJECT>_dyn_mod.mat
%
% AUTHOR: Mohammad Reza Salehi
% =========================================================================

clc;
clear;
close all;

%% 1. CONFIGURATION
% =========================================================================

project_root = 'F:\PhD Code\My_PhD_Project';

% -------------------------------------------------------------------------
% DATA TYPE
% -------------------------------------------------------------------------
% Options:
%   'inter_parcellated'
%   'processed_clean'

further_cleaning = 'inter_parcellated';

% -------------------------------------------------------------------------
% PIPELINE NAME
% -------------------------------------------------------------------------
% Only relevant for processed_clean

cleaning_GSR = 'With_GSR';

% -------------------------------------------------------------------------
% SESSIONS
% -------------------------------------------------------------------------

sessions = { ...
    'REST1_LR', ...
    'REST1_RL', ...
    'REST2_LR', ...
    'REST2_RL'};

% -------------------------------------------------------------------------
% BASE RESULT DIRECTORY
% -------------------------------------------------------------------------

base_result_dir = fullfile(project_root, ...
    'results', ...
    'Modulated_FC', ...
    'retest');

%% 2. MAIN LOOP
% =========================================================================

fprintf('======================================================\n');
fprintf('RUNNING RETEST DYNAMIC MODULATION ANALYSIS\n');
fprintf('======================================================\n');

t_all = tic;

for s = 1:length(sessions)

    session = sessions{s};

    fprintf('\n==================================================\n');
    fprintf('SESSION: %s\n', session);
    fprintf('==================================================\n');

    %% Define Input/Output Paths
    % =====================================================================

    if strcmp(further_cleaning, 'processed_clean')

        data_root = fullfile(project_root, ...
            'data', ...
            'processed_clean', ...
            'Schaefer200_Kong17', ...
            'retest', ...
            session, ...
            cleaning_GSR);

        output_dir = fullfile(base_result_dir, ...
            session, ...
            cleaning_GSR);

    elseif strcmp(further_cleaning, 'inter_parcellated')

        data_root = fullfile(project_root, ...
            'data', ...
            'inter_parcellated', ...
            'Schaefer200_Kong17', ...
            'retest', ...
            session);

        output_dir = fullfile(base_result_dir, ...
            session, ...
            'No_clean');

    else
        error('Unknown further_cleaning option.');
    end

    % Create output directory
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    %% Subject List
    % =====================================================================

    subjects = dir(fullfile(data_root, '*.mat'));

    fprintf('Found %d files.\n', length(subjects));

    %% Waitbar
    % =====================================================================

    h = waitbar(0, ...
        sprintf('Session %s Initialization...', session));

    %% Subject Loop
    % =====================================================================

    for i = 1:length(subjects)

        subj_name = subjects(i).name;

        subj_path = fullfile(subjects(i).folder, subj_name);

        fprintf('[%s] Processing %s\n', session, subj_name);

        % -----------------------------------------------------------------
        % Load Data
        % -----------------------------------------------------------------

        tmp = load(subj_path);

        if strcmp(further_cleaning, 'processed_clean')

            ts = tmp.clean_data';

        elseif strcmp(further_cleaning, 'inter_parcellated')

            ts = tmp.data_on_atlas';

        end

        % -----------------------------------------------------------------
        % CORE CALCULATION
        % -----------------------------------------------------------------

        subj_metrics = calc_dynamic_modulation_metrics( ...
            ts, ...
            'n_levels', 20, ...
            'req_metrics', {'all'});

        % -----------------------------------------------------------------
        % Save Output
        % -----------------------------------------------------------------

        [~, name_no_ext, ~] = fileparts(subj_name);

        save_path = fullfile(output_dir, ...
            [name_no_ext '_dyn_mod.mat']);

        save(save_path, '-struct', 'subj_metrics', '-v7.3');

        % -----------------------------------------------------------------
        % Update Waitbar
        % -----------------------------------------------------------------

        waitbar(i / length(subjects), ...
            h, ...
            sprintf('Session %s | %d/%d | %s', ...
            session, i, length(subjects), name_no_ext));

    end

    close(h);

    fprintf('\nSession %s completed.\n', session);
    fprintf('Results saved to:\n%s\n', output_dir);

end

%% 3. FINISH
% =========================================================================

fprintf('\n======================================================\n');
fprintf('ALL RETEST ANALYSES COMPLETED\n');
fprintf('Total elapsed time: %.2f minutes\n', toc(t_all)/60);
fprintf('======================================================\n');