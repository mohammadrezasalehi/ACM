function [roi_names, roi_nets] = load_atlas_info(atlas_info_path)
% LOAD_ATLAS_INFO  Load ROI names and network labels from atlas info file
%
%   [roi_names, roi_nets] = load_atlas_info(atlas_info_path)
%
%   Inputs:
%       atlas_info_path : path to atlas info text file
%
%   Outputs:
%       roi_names : cell array of ROI/parcel names
%       roi_nets  : cell array of corresponding network labels

    fid = fopen(atlas_info_path, 'r');
    if fid == -1
        error('Cannot open atlas info file: %s', atlas_info_path);
    end

    roi_names = {};
    roi_nets  = {};
    idx = 1;

    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line)
            continue;
        end

        if isletter(line(1)) || contains(line, '17networks')
            roi_names{idx} = line;

            parts = strsplit(line, '_');
            if numel(parts) >= 3
                roi_nets{idx} = parts{3};
            else
                roi_nets{idx} = 'Unk';
            end

            % skip next line (as in original code)
            fgetl(fid);
            idx = idx + 1;
        end
    end

    fclose(fid);
end
