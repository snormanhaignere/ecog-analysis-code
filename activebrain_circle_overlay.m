function activebrain_circle_overlay(exp, subjid, hemi, electrode_activations, electrode_indices, varargin)

% requires activeBrain and export_fig_v2

% default parameters
I.radii = 2.5;
I.figname = '';
I.cmap = [];
I.cmap_range = [min(electrode_activations), max(electrode_activations)];
I.fading = 0;
I.only_plot_given_electrodes = false;
I.figh = matlab.ui.Figure.empty;
I = parse_optInputs_keyvalue(varargin, I);

% load brain model in talaraich coordinates
load([root_directory '/' exp '/data/Anatomical/' subjid '/BrainModel.mat'], ...
    'cortex', 'viewstruct', 'cmapstruct', 'ix', 'tala', 'vcontribs');

% viewing parameters
viewstruct.enablelight = 1;
viewstruct.what2view={'brain','trielectrodes'};
switch hemi
    case 'rh'
        viewstruct.lightpos = [150, 0, 0]; % light
        viewstruct.viewvect = [90, 0]; % view
    case 'lh'
        viewstruct.lightpos = [-150, 0, 0];
        viewstruct.viewvect = [270, 0]; % view
    otherwise
        error('hemi should be rh or lh not %s', hemi);
end

% if trielectrodes are not present, create them
if ~isfield(tala, 'trielectrodes') %#ok<NODEF>
    tala.trielectrodes = tala.electrodes;
end

% set activations
if I.only_plot_given_electrodes
    tala.activations = electrode_activations';
    tala.electrodes = tala.electrodes(electrode_indices,:);
    tala.trielectrodes = tala.trielectrodes(electrode_indices,:);
    vcontribs_new = [];
    for i = 1:length(vcontribs) %#ok<NODEF>
        elec = vcontribs(i).contribs(:,2);
        [~,ai,bi] = intersect(elec, electrode_indices);
        if ~isempty(ai)
            vcontribs_new = [vcontribs_new, vcontribs(i)]; %#ok<AGROW>
            vcontribs_new(end).contribs = vcontribs(i).contribs(ai,:);
            vcontribs_new(end).contribs(:,2) = bi; 
        end
    end
    vcontribs = vcontribs_new;
else
    tala.activations = zeros(size(tala.electrodes,1),1);
    [~,xi] = intersect(electrode_indices, 1:size(tala.electrodes,1));
    tala.activations(electrode_indices(xi)) = electrode_activations(xi);
    clear xi;
end

% sets min and max of colorbar
cmapstruct.cmin = I.cmap_range(1);
cmapstruct.cmax = I.cmap_range(2);
cmapstruct.enablecolorbar = 1;
cmapstruct.fading = I.fading;

% color map to use
if ~isempty(I.cmap)
    cmapstruct.cmap = I.cmap;
end

% plot
if isempty(I.figh)
    I.figh = figure;
end
clf(I.figh);
set(I.figh, 'Position', [0 0 1000 600]);
activateBrain_colorelectrodes(cortex, vcontribs, tala, ix, ...
    cmapstruct, viewstruct, tala.activations, I.radii);

% save
if ~isempty(I.figname)
    box off;    
    set(gcf, 'PaperSize', [10 6]);
    set(gcf, 'PaperPosition', [0.25 0.25 9.5 5.5]);
    print([I.figname '.pdf'],'-dpdf');
    print([I.figname '.png'],'-dpng', '-r200');
end