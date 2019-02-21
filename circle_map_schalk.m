function [map, elec_found] = circle_map_schalk(subjid, elec_inds, elec_vals, varargin)

% 2019-01-30: Creates a map with electrode values represented as circles,
% used for schalk data

global root_directory;
freesurfer_directory = [root_directory '/freesurfer'];
addpath(genpath([root_directory '/fmri-analysis-beta/code']));

n_verts_fsaverage = 163842;

% radius in mm
I.radmm = 2;
I.map = nan(n_verts_fsaverage, 2);
I.locstring = '';
I = parse_optInputs_keyvalue(varargin, I);

assert(length(elec_inds)==length(elec_vals));
n_elec = length(elec_vals);
hemis = {'rh', 'lh'};
elec_found = false(1,n_elec);
map = I.map;
for q = 1:length(hemis)
    
    hemi = hemis{q};
    
    for j = 1:n_elec
        
        % find mean vertex as weighted average of inflated vertex locations
        fsaverage_elec_directory = [freesurfer_directory '/myfsaverage/elec-location-maps' I.locstring '/' subjid];
        map_fsaverage = [fsaverage_elec_directory '/' hemi '.elec' num2str(elec_inds(j)) '.mgh'];
        br = MRIread(map_fsaverage);
        w = br.vol(:);
        [~, wi] = max(w(:));
        [inflated_verts, ~] = freesurfer_read_surf([freesurfer_directory '/myfsaverage/surf/' hemi '.inflated']);
        max_vert = inflated_verts(wi,:);
        
        % create circle
        d = sqrt(sum(bsxfun(@minus, inflated_verts, max_vert).^2, 2));
        map(d<I.radmm, q) = elec_vals(j);
        clear yi;
        
    end
end

