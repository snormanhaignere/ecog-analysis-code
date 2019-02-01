function E = read_elec_coords_pial(subjid, hemi)

% Reads the electrode coordinate information on the pial surface
% 
% 2019-01-18: Created, Sam NH

global root_directory;

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