function [STRFs, params] = calcSTRFs_out(out,channels,gap,compressresp,plotvals,useSTRFLab,nC,nDS_freq,nDS_time,lags)
% [STRFs] = calcSTRFs_out(out,channels,gap,compressresp,plotvals,useSTRFLab,nC,nDS_freq,nDS_time)
%
% Wrapper for STRF- both James?s and STRFLab using STRFCrossValidate
% First extracts data (in cell for STRFlab, in array for James), then runs STRFs for each one
% Output is struct that is roughly aligned with each other
% Will plot r_values if desired
% Can be STRF lab or not (STRFlab is the default, having weird issues with Nima?s code
%
% Written by Laura, NAPLab March 2017

if ~exist('channels','var') || isempty(channels)
    channels = 1:size(out(1).resp,1);
end
if ~exist('gap','var') || isempty(gap)
    gap = 0;
end
if ~exist('plotvals','var') || isempty(plotvals)
    plotvals = 0;
end
if ~exist('useSTRFLab','var') || isempty(useSTRFLab)
    useSTRFLab = 1;
    dispDefaultMessage(useSTRFLab,'useSTRFLab');
end
if ~exist('nC','var') || isempty(nC)
    nC = 0.33;
    compressstim = 1; 
    dispDefaultMessage(nC,'nC');
elseif nC == 1
    compressstim = 0;
else
    compressstim = 1;
end
if ~exist('nDS_freq','var') || isempty(nDS_freq)
    nDS_freq = 4; 
    resamp_freq = 1;
    dispDefaultMessage(nDS_freq,'nDS_freq');
elseif nC == 1
    resamp_freq = 0;
else
    resamp_freq = 1;
end
if ~exist('nDS_time','var') || isempty(nDS_time)
    nDS_time = 1;
    resamp_time = 0;
elseif nDS_time == 1
    resamp_time = 0;
else
    resamp_time = 1;
end
if ~exist('compressresp','var') || isempty(compressresp)
    compressresp = 1; 
    dispDefaultMessage(compressresp,'compressresp');
end
if ~exist('lags','var') || isempty(lags)
    lags = [];
end

fsdata = out(1).dataf;
% Check whether stim is mono
removestimchan = size(out(1).aud,3) ~= 1;
if removestimchan
    warning('stim is not mono; keeping one channel only');
    whichstimchan = input('Which stim chan should be kept? ');
end


if useSTRFLab
    
    disp('...extracting stim and resp');
    
    % Decide about splitting up the trials if there are too few for the crossvalidations
    numstim = length(out);
    splitstim = numstim < 200;
    if splitstim
        numsplit = input('How many segments should each trial be split into? ');
    else
        numsplit = 1;
    end
    
    % Set up variables
    resp = cell(1,numstim*numsplit); stim = cell(1,numstim*numsplit);
    
    ind = 1;
    for i = 1:numstim
        
        % Grab stim and resp
        thisresp = mean(out(i).resp,3);
        if removestimchan
            thisstim = out(i).aud(:,:,whichstimchan);
        else
            thisstim = out(i).aud;
        end
        
        % Compressions if applicable
        if compressresp
            thisresp = mapstd(thisresp);
            thisresp = 10*tanh(thisresp/10);
        end
        if compressstim
            thisstim = thisstim.^nC;
        end
        
        % Resampling if applicable
        if resamp_time
            thisstim = nt_dsample(thisstim,nDS_time);
            thisresp = nt_dsample(thisresp,nDS_time);
        end
        if resamp_freq
            thisstim = nt_dsample(thisstim,nDS_freq);
        end
        
        
        % Concatenate responses into cell arrays
        if numsplit == 1
            resp{ind} = thisresp;
            stim{ind} = thisstim;
            ind = ind+1;
        else
            windowsize = floor(size(thisresp,2)/numsplit);
            for j = 1:numsplit
                if j == numsplit
                    thiswindow = (j-1)*windowsize:size(thisresp,2);
                else
                    thiswindow = (j-1)*windowsize+1:j*windowsize;
                end
                resp{ind} = thisresp(:,thiswindow);
                stim{ind} = thisstim(:,thiswindow);
                ind = ind+1;
            end
        end
        
        
    end
    
    
    % Run STRF code
    disp('...running STRFCrossValidate');
    defaultparams = input('Use default tolerance and sparsity parameters? ');
    if defaultparams
        tolval=[.005 .01 .05 .1] ; sparseval=[4 6 8 16 32];
    else
        tolval = input('Input desired tolerance values: ');
        sparseval = input('Input desired sparsity values: ');
    end
    for i = 1:length(channels)
        [strf, params, presp, nresp, corrs] = STRFCrossValidate(stim,resp,channels(i),40,fsdata,{0,{length(lags),lags},'DirectFit'},{sparseval,tolval});
        STRFs.strf(:,:,:,i) = permute(strf,[3 1 2]); % should look at rearranging this to be more like what James uses
        STRFs.r_values(:,i) = corrs; disp(mean(corrs));
        STRFs.prediction(i,:) = presp;
        STRFs.actual(i,:) = nresp;
        STRFs.params(i) = params;
    end
    
    if plotvals
        figure; plot(mean(STRFs.r_values),'o-'); xlabel('Electrode'); ylabel('R Value');
    end
    
else
    
    disp('...extracting stim and resp');
    
    % Set up variables
    resp = []; stim = [];
    
    % Loop to load stim and resp into
    for i = 1:length(out)
        if removestimchan
            thisstim = out(i).aud(:,:,whichstimchan);
        else
            thisstim = out(i).aud;
        end
        thisresp = mean(out(i).resp,3);
        stim = cat(2,stim,thisstim);
        resp = cat(2,resp,thisresp);
        stim = cat(2,stim,zeros(size(stim,1),gap));
        resp = cat(2,resp,zeros(size(resp,1),gap));
    end
    
    % Compress resp and stim
    if compressresp
        resp = mapstd(resp);
        resp = 10*tanh(resp/10);
    end
    
    % Settings
    disp('...running Run_STRF');
    k = 10;
    method = 'ridge';
    
    resp = resp(channels,:);
    STRFs = Run_STRF(stim,resp,1,lags,nDS_freq,nDS_time,nC,k,method);
    if plotvals
        figure; plot(mean(STRFs.r_values),'o-'); xlabel('Electrode'); ylabel('R Value');
    end
end

params.channels = channels;
params.gap = gap;
params.compressresp = compressresp;
params.plotvals = plotvals;
params.useSTRFLab = useSTRFLab;
params.nC = nC;
params.nDS_freq = nDS_freq;
params.nDS_time = nDS_time;
if useSTRFLab
    params.tolval = tolval;
    params.sparseval = sparseval;
end

end