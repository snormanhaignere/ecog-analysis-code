function [map, elec_found] = circle_map(subjid, elec_vals, varargin)

% 2019-01-30: Creates a map with electrode values represented as circles

global root_directory;
freesurfer_directory = [root_directory '/freesurfer'];
addpath(genpath([root_directory '/fmri-analysis-beta/code']));

n_verts_fsaverage = 163842;

% radius in mm
I.radmm = 2;
I.map = nan(n_verts_fsaverage, 2);
I.elec_names = read_channel_names(subjid);
I.elec_inds = 1:length(I.elec_names);
I.locstring = '';
[I, C] = parse_optInputs_keyvalue(varargin, I);

if C.elec_inds && ~C.elec_names
    I.elec_names = I.elec_names(I.elec_inds);
end

if C.elec_names && ~C.elec_inds
    all_channel_names = read_channel_names(subjid);
    I.elec_inds = nan(size(I.elec_names));
    for i = 1:length(I.elec_names)
        I.elec_inds(i) = find(ismember(all_channel_names, I.elec_names{i}));
    end
end

assert(length(I.elec_names)==length(elec_vals));
n_elec = length(elec_vals);
hemis = {'rh', 'lh'};
elec_found = false(1,n_elec);
map = I.map;
for q = 1:length(hemis)
    
    hemi = hemis{q};
    
    % load electrode info
    E = read_elec_coords_pial(subjid, hemi);
    
    if isempty(E)
        continue;
    end
    
    for j = 1:length(E.chnames)
        
        xi = find(ismember(I.elec_names, E.chnames{j}));
        
        if ~isempty(xi)
                        
            % subject-specific map
            assert(length(xi)==1);
            
            % indicate that the electrode has been found
            elec_found(xi) = true;
            
            % find mean vertex as weighted average of inflated vertex locations
            fsaverage_elec_directory = [freesurfer_directory '/myfsaverage/elec-location-maps' I.locstring '/' subjid];
            map_fsaverage = [fsaverage_elec_directory '/' hemi '.elec' num2str(I.elec_inds(xi)) '-' I.elec_names{xi} '.mgh'];
            br = MRIread(map_fsaverage);
            w = br.vol(:);
            [~, wi] = max(w(:));
            [inflated_verts, ~] = freesurfer_read_surf([freesurfer_directory '/myfsaverage/surf/' hemi '.inflated']);
            max_vert = inflated_verts(wi,:);
            
            % create circle
            d = sqrt(sum(bsxfun(@minus, inflated_verts, max_vert).^2, 2));
            map(d<I.radmm, q) = elec_vals(xi);
            clear yi;
            
        end
        clear xi;
    end
end

