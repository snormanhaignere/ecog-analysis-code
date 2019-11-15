function [ reconout ] = createReconOut(recon, out, reconparams, tag)
% [ reconout ] = createReconOut( recon, out )
% Turns recon struct into an out structure (so it can be processed in
% units)

if ~exist('tag','var') || isempty('tag')
    tag = '';
end


%% Check that recon and out are the same length
[stim, ~] = extractRespStim_out(out);

if ~isequal(recon.stim_test,stim)
    [stim2, ~] = extractRespStim_out(out,[],'cell');
    if ~isequal(recon.stim_test,stim2)
        error('ERROR: recon and out stimuli do not match.');
    end
end

%% Deal rstim into out structure

reconout = out;
reconparams.g = recon.g;
reconparams.corr = recon.corr;

start = 1;
for i = 1:length(out)
    thislength = size(out(i).aud,2)-1;
    if iscell(recon.rstim)
        reconout(i).(['rstim' tag]) = recon.rstim{i};
    else
        reconout(i).(['rstim' tag]) = recon.rstim(:,start:start+thislength);
    end
    reconout(i).(['reconparams' tag]) = reconparams;
    start = start+thislength+1;
end

if start ~= size(recon.rstim,2)+1
    warning('Warning: Not all of stim was used!');
end


end

