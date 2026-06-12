%% Calculates centroids by reading GIFTI files directly in MATLAB.
clc; clear;

%% 1. CONFIGURATION
% addresses
project_root = 'F:\PhD Code\My_PhD_Project';
% add toolbox to read cifti and gifti files
addpath(genpath(fullfile(project_root, ...
    'tools', 'cifti-matlab-master')));
% Workbench
wb_cmd_path = fullfile(project_root, ...
    'tools', 'workbench', 'bin_windows64', 'wb_command.exe');
% Atlas
dlabel_path = fullfile(project_root, 'data', 'atlases', ...
    'Schaefer2018_200Parcels_Kong2022_17Networks_order.dlabel.nii');
% Surfaces
surf_L_path = fullfile(project_root, 'data', 'atlases', 'geometry_data', ...
    'S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii');
surf_R_path = fullfile(project_root, 'data', 'atlases', 'geometry_data', ...
    'S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii');
% Output
output_path = fullfile(project_root, 'data', 'atlases', 'centroids',...
    'centroids_schaefer200_kong.mat');

%% -------- 2. CHECK FILES --------

assert(exist(wb_cmd_path,'file')==2, 'wb_command not found');
assert(exist(dlabel_path,'file')==2, 'dlabel file not found');
assert(exist(surf_L_path,'file')==2, 'Left surface not found');
assert(exist(surf_R_path,'file')==2, 'Right surface not found');

%% 2. SEPARATE HEMISPHERES

label_L_gii = 'temp_L.label.gii';
label_R_gii = 'temp_R.label.gii';

fprintf('Splitting CIFTI to GIFTI...\n');
system(sprintf('"%s" -cifti-separate "%s" COLUMN -label CORTEX_LEFT "%s" -label CORTEX_RIGHT "%s"', ...
    wb_cmd_path, dlabel_path, label_L_gii, label_R_gii));

%% 3. LOAD DATA IN MATLAB (USING GIFTI LIBRARY)
fprintf('Loading GIFTI surfaces and labels...\n');

try
    % Load geometry (XYZ coordinates)
    S_L = gifti(surf_L_path); 
    S_R = gifti(surf_R_path);
    
    % Load labels (which vertex belongs to which region)
    L_L = gifti(label_L_gii);
    L_R = gifti(label_R_gii);
    
    % Extract raw data
    verts_L = S_L.vertices; % Matrix [32492 x 3]
    verts_R = S_R.vertices; % Matrix [32492 x 3]
    
    labels_L = L_L.cdata;   % Vector [32492 x 1]
    labels_R = L_R.cdata;   % Vector [32492 x 1]
    
catch ME
    error('Gifti library not found or file error. Ensure FieldTrip/SPM is in path.');
end

%% 4. CALCULATE CENTROIDS
fprintf('Calculating Centroids...\n');

n_rois = 200;
centroids = zeros(n_rois, 3);

for i = 1:n_rois
    % Schiffer Atlas: 1-100 left, 101-200 right
    if i <= 100
        % Find vertices belonging to region i in the left hemisphere
        idx = find(labels_L == i);
        if ~isempty(idx)
            % Average coordinates of these vertices
            centroids(i, :) = mean(verts_L(idx, :), 1);
        else
            centroids(i, :) = [NaN NaN NaN];
        end
    else
        % Right hemisphere
        idx = find(labels_R == i);
        if ~isempty(idx)
            centroids(i, :) = mean(verts_R(idx, :), 1);
        else
            centroids(i, :) = [NaN NaN NaN];
        end
    end
end

%% 5. SAVE & CLEANUP
save(output_path, 'centroids');
fprintf('Done! Centroids saved to centroids_schaefer200_kong.mat\n');

% Clear temporary files
delete('temp_L.label.gii');
delete('temp_R.label.gii');

% Show a few examples to be sure
fprintf('\nCheck Results:\n');
fprintf('Node 1 (Left):  [%f, %f, %f]\n', centroids(1,:));
fprintf('Node 101 (Right): [%f, %f, %f]\n', centroids(101,:));