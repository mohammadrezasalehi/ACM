%% =========================================================================
% SCRIPT: find_retest_subject_indices.m
% =========================================================================
% DESCRIPTION:
%   Finds which indices the 8 retest subjects occupy inside the
%   original 100 unrelated subject list.
%
% AUTHOR:
%   Mohammad Reza Salehi
% =========================================================================

clc;
clear;
close all;

%% 1. PATHS
% =========================================================================

retest_dir = ...
    'F:\PhD Code\My_PhD_Project\data\inter_parcellated\Schaefer200_Kong17\retest\REST1_LR';

main_dir = ...
    'F:\PhD Code\My_PhD_Project\data\inter_parcellated\Schaefer200_Kong17\REST1_LR';

%% 2. GET FILE LISTS
% =========================================================================

% Retest subjects
retest_files = dir(fullfile(retest_dir, '*.mat'));

% Main 100 unrelated subjects
main_files = dir(fullfile(main_dir, '*.mat'));

%% 3. EXTRACT SUBJECT IDs
% =========================================================================

retest_ids = cell(length(retest_files),1);

for i = 1:length(retest_files)

    [~, name, ~] = fileparts(retest_files(i).name);

    % Example:
    % 103818.mat -> 103818
    retest_ids{i} = name;

end

main_ids = cell(length(main_files),1);

for i = 1:length(main_files)

    [~, name, ~] = fileparts(main_files(i).name);

    main_ids{i} = name;

end

%% 4. FIND MATCHING INDICES
% =========================================================================

fprintf('\n====================================================\n');
fprintf('Matching Retest Subjects Inside Main Dataset\n');
fprintf('====================================================\n\n');

matched_indices = zeros(length(retest_ids),1);

for i = 1:length(retest_ids)

    idx = find(strcmp(main_ids, retest_ids{i}));

    if isempty(idx)

        fprintf('Subject %s --> NOT FOUND\n', retest_ids{i});
        matched_indices(i) = NaN;

    else

        fprintf('Subject %s --> Index %d\n', ...
            retest_ids{i}, idx);

        matched_indices(i) = idx;

    end
end

%% 5. SUMMARY
% =========================================================================

fprintf('\n====================================================\n');
fprintf('Matched Indices Vector:\n');
fprintf('====================================================\n');

disp(matched_indices');

%% Optional save
save('retest_subject_indices.mat', ...
    'matched_indices', ...
    'retest_ids');