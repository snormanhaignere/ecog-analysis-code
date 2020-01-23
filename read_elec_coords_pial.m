function E = read_elec_coords_pial(subjid, hemi)

% Reads the electrode coordinate information on the pial surface
%
% 2019-01-18: Created, Sam NH

global root_directory;

% check if coordinate file created by iELVis is present
% if no assume the coordinate files were created using ntools_elec
% (probably provided by Hugh)
clear E;
coord_file = [root_directory '/freesurfer/' subjid '/elec_recon/POSTIMPLANT.txt'];
if exist(coord_file, 'file')
    
    fid = fopen(coord_file, 'r');
    x = textscan(fid, '%f%f%f', 'HeaderLines', 2); fclose(fid);
    E.coords = cat(2,x{:});
    clear coord_file fid x;
    
    chname_file = [root_directory '/freesurfer/' subjid '/elec_recon/electrodeNames.txt'];
    fid = fopen(chname_file, 'r');
    x = textscan(fid, '%s%s%s', 'HeaderLines', 2); fclose(fid);
    E.chnames = x{1};
    E.elec_type = x{2};
    E.hemi = x{3};
    clear x chname_File fid;
        
    % change format of hemi
    E.hemi = strrep(E.hemi, 'R', 'rh');
    E.hemi = strrep(E.hemi, 'L', 'lh');
    
    switch hemi
        case {'rh', 'lh'}
            xi = ismember(E.hemi, hemi);
            E.coords = E.coords(xi, :);
            E.chnames = E.chnames(xi);
            E.elec_type = E.elec_type(xi);
            E.hemi = E.hemi(xi);
            clear xi;
        case 'both'
            % do nothing
        otherwise
            error('%s doesn''t match format\n', hemi);
    end
    
else
    
    % read electrode information
    elec_directory = [root_directory '/freesurfer/' subjid '/elec'];
    if strcmp(hemi, 'both')
        pial_elec_coords_file = mydir(elec_directory, 'coor_T1_');
        file_found = false;
        for i = 1:length(pial_elec_coords_file)
            if isempty(strfind(pial_elec_coords_file{i}, '_LH_')) ...
                    && isempty(strfind(pial_elec_coords_file{i}, '_RH_'))
                pial_elec_coords_file = pial_elec_coords_file(i);
                file_found = true;
                break;
            end
        end
        if ~file_found
            error('File not found');
        end
    else
        pial_elec_coords_file = mydir(elec_directory, ['coor_T1_' upper(hemi)]);
    end
    
    % if file not present return nothing
    if isempty(pial_elec_coords_file)
        E = [];
        return;
    end
    
    % otherwise return the info
    assert(length(pial_elec_coords_file)==1);
    fid = fopen([elec_directory '/' pial_elec_coords_file{1}]);
    x = textscan(fid, '%s%f%f%f%s'); clear fid;
    E.chnames = x{1};
    E.coords = cat(2, x{2:4});
    E.elec_type = x{5};
    clear x;
    
end

% for some subjects append hemi to channel name
switch subjid
    case {'084_CUBF42', '085_CU102', '083_CU101'}
        for i = 1:length(E.chnames)
            if strcmp(subjid, '085_CU102') && strcmp(E.chnames{i}(1:4), 'LAMY')   
                % do nothing
            else
                switch E.hemi{i}
                    case 'rh'
                        E.chnames{i} = ['R' E.chnames{i}];
                    case 'lh'
                        E.chnames{i} = ['L' E.chnames{i}];
                    otherwise
                        error('No matching hemi');
                end
            end
        end
end