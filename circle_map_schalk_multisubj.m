function [map, hemi_assignment] = circle_map_schalk_multisubj(all_subjid, ei, si, evals, varargin)

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
I.order = 'ascend';
I.keyboard = false;
I = parse_optInputs_keyvalue(varargin, I);

% debug
if I.keyboard
    keyboard;
end

switch I.order
    case 'ascend'
        [~, xi] = sort(evals);
    otherwise
        error('Switch fell through');
end
ei = ei(xi);
si = si(xi);
evals = evals(xi);
clear xi;
assert(length(ei)==length(evals));
assert(length(ei)==length(si));
n_elec = length(evals);

% inflated vertices on the fsaverage brain
hemis = {'rh', 'lh'};
inflated_verts = cell(1,2);
for q = 1:2
    [inflated_verts{q}, ~] = freesurfer_read_surf([freesurfer_directory '/myfsaverage/surf/' hemis{q} '.inflated']);
end
inflated_verts_both_hemi = [inflated_verts{1}; inflated_verts{2}];
hemi_labels = [ones(size(inflated_verts{1},1),1); 2*ones(size(inflated_verts{2},1),1)];

map = I.map(:);
hemi_assignment = nan(1, n_elec);
for j = 1:n_elec

    % find mean vertex as weighted average of inflated vertex locations
    fsaverage_elec_directory = [freesurfer_directory '/myfsaverage/elec-location-maps-constrained/' I.locstring '/' all_subjid{si(j)}];
    w = [];
    for q = 1:2
        map_fsaverage = [fsaverage_elec_directory '/' hemis{q} '.elec' num2str(ei(j)) '.mgh'];
        br = MRIread(map_fsaverage);
        w = [w; br.vol(:)]; %#ok<*AGROW>
    end
    [~, wi] = max(w(:));
    max_vert = inflated_verts_both_hemi(wi,:);
    hemi_assignment(j) = hemi_labels(wi);
    
    % create circle
    d = sqrt(sum(bsxfun(@minus, inflated_verts_both_hemi, max_vert).^2, 2));
    map(d<I.radmm) = evals(j);
        
end

% separate out dimensions
map = reshape(map, [n_verts_fsaverage, 2]);

