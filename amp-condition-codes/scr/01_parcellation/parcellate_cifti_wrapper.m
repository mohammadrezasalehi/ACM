function [success, msg] = parcellate_cifti_wrapper(input_nii_path, atlas_path, output_mat_path, wb_cmd_path)
% =========================================================================
% FUNCTION: parcellate_cifti_wrapper
% =========================================================================
% DESCRIPTION:
%   Wraps the Connectome Workbench command to parcellate a CIFTI file and
%   converts the result immediately to a .mat file to save space.
%
% INPUTS:
%   input_nii_path  : Full path to the input .dtseries.nii file
%   atlas_path      : Full path to the .dlabel.nii atlas file
%   output_mat_path : Full path where the .mat file should be saved
%   wb_cmd_path     : Path to wb_command executable
%
% OUTPUTS:
%   success : Boolean (true if successful, false otherwise)
%   msg     : Error message or status report
% =========================================================================

    success = false;
    msg = '';

    % 1. Validation
    if ~exist(input_nii_path, 'file')
        msg = 'Input file does not exist.';
        return;
    end
    
    % Temporary file path (created in the same folder as output, usually)
    [out_dir, out_name, ~] = fileparts(output_mat_path);
    temp_ptseries = fullfile(out_dir, [out_name, '_temp.ptseries.nii']);

    try
        % 2. Execute wb_command
        % Syntax: -cifti-parcellate <cifti-in> <cifti-label> COLUMN <cifti-out>
        cmd = sprintf('"%s" -cifti-parcellate "%s" "%s" COLUMN "%s"', ...
            wb_cmd_path, input_nii_path, atlas_path, temp_ptseries);
        
        [status, cmdout] = system(cmd);
        
        if status ~= 0
            msg = sprintf('wb_command failed: %s', cmdout);
            % Clean up if partial file created
            if exist(temp_ptseries, 'file'), delete(temp_ptseries); end
            return;
        end
        
        % 3. Convert to .mat (using cifti-matlab/FieldTrip)
        cifti_data = ft_read_cifti(temp_ptseries);
        
        % Check fields to extract data matrix safely
        if isfield(cifti_data, 'ptseries')
            data_on_atlas = cifti_data.ptseries;
        elseif isfield(cifti_data, 'cdata')
            data_on_atlas = cifti_data.cdata;
        else
            data_on_atlas = cifti_data.dtseries; % Fallback
        end
        
        % 4. Save
        save(output_mat_path, 'data_on_atlas');
        
        % 5. Cleanup
        delete(temp_ptseries);
        
        success = true;
        msg = 'Processed successfully';
        
    catch ME
        msg = sprintf('MATLAB Error: %s', ME.message);
        if exist(temp_ptseries, 'file'), delete(temp_ptseries); end
    end
end