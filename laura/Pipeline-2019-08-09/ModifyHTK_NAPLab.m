%% User Inputs: Set parameters, choose what actions to take
% Standardized for Lab Use by Laura Long 06/2019

% Enter subject ID.
subject = '086_NY723'; % subject name
task = 'Scrambling'; % task name

% Enter dpath. Should be the patient folder, containing 'original' directory.
% If empty, defaults to current directory.
dpath = ['/Users/LauraLong/Documents/Lab/ECoGData/' subject];

% Enter recording system. Necessary for ConvertToHTK, otherwise optional.
% Options: 'tdt' 'blackrock' 'edf' 'eeg' 'newtdt'
recordingsystem = 'blackrock';

% Enter which data you want to start with.
% If empty, defaults to htkraw.
whichdata = '';

% Enter additional parameters.
% If empty, defaults to all blocks, all channels, and 1000 Hz.
blocks = 1; % blocks to process
channelinds = []; % specify which electrodes to include (if eliminating some)
dataf = []; % desired fs to save data

% Possible data processing steps. Set to 1 to take each action.
% Will execute in order listed.
converttohtk = 1; % converts data from original format to .htk
writechnames = 0; % manually writes channel names (if they aren't already available in original file; rare)
findcleanvars = 0; % runs data checks for bad channels, referencing, filtering, etc; generates plots and allows user to pick best params
cleandata = 0; % applies data cleaning found from findcleanvars (can also input params directly without running findcleanvars)
highgamma = 0; % extracts and saves high gamma 
freqband = 0; % extracts and saves other frequency bands
makestim = 1; % finds and writes StimOrder to file; best to run one block at a time
findevnt_aud = 0; % finds stimulus alignment where audio was recorded
findevnt_trig = 1; % finds stimulus alignment where triggers were recorded
genout = 0; % generates out structure 

% Extra params are necessary for several functions; check below! 
plotresults = 1; % if 1/true, will generate several reporting plots along the way


%%  Convert to HTK

if converttohtk % converts data from original format to .htk
    ConvertToHTK_NAPLab(dpath,recordingsystem,blocks,channelinds,dataf);
end

if writechnames % manually writes channel names; only used if names are not available (or incorrect) in original format 
    trodenames = {}; % channel name prefixes
    trodelens = []; % number of channels for each prefix
    WriteChNames_NAPLab(dpath,blocks,trodenames,trodelens);
end


%% Clean data

if findcleanvars
    cleanblock = 1;
    cleanvars = FindCleanDataVariables_NAPLab(dpath,cleanblock,task);
    v2struct(cleanvars);
end

if cleandata
    
    
    % If you didn't just find the cleanvars, ask about loading from file
    if ~findcleanvars
        dataqualityfile = dir([dpath filesep 'processed' filesep 'B' num2str(cleanblock) filesep  '*dataqualityvars.mat']);
        dataqualityfile = [dpath filesep 'processed' filesep 'B' num2str(cleanblock) filesep dataqualityfile.name];
        if exist(dataqualityfile,'file');
            loadinput = input('Load cleanvars from file for HTKtoCleanData? ');
            if loadinput
                disp('...loading cleanvars from file');
                load(dataqualityfile);
                v2struct(cleanvars);
            end
        end
    end
    
    % Clean Data!
    for i = 1:length(blocks)
        HTKtoCleanData_NAPLab(dpath,whichdata,blocks(i),channelinds,dataf,filters.highpassfreq,reftype,amplifiers,channels.bad,epochs,channels.iffy,normdata,plotresults);
    end
    
    % Update whichdata flag
    if isempty(whichdata)
        whichdata = 'cleaned';
    else
        whichdata = [whichdata '_cleaned'];
    end
    
end

%% Extract Frequency Bands

if highgamma
    
    HTKtoHighGamma_NAPLab(dpath,whichdata,blocks,channelinds,dataf);
    
    % Update whichdata flag
    if isempty(whichdata)
        whichdata = 'highgamma';
    else
        whichdata = [whichdata '_highgamma'];
    end
end

if freqband
    freqRange = [1 70];
    getenvelope = 0;
    datatag = 'lowfreq';
    
    HTKtoFreqBand_NAPLab(dpath,whichdata,blocks,channelinds,dataf,freqRange,getenvelope,datatag);
    
    if isempty(whichdata)
        whichdata = datatag;
    else
        whichdata = [whichdata '_' datatag];
    end
end

%% Make Stimulus Folders

if makestim
    
    % Inputs for this function only
    stimpath = '/Users/LauraLong/Documents/Lab/ECoGData/Task_Sounds';
    trials = []; % for fixed order, enter array indicating which trial; if fixed order, leave empty and code will look for StimOrder in behavioral/B#
    portion = {}; % 
    
    % Generate stim folders!
    makeStimFolder_NAPLab(dpath,blocks,stimpath,task,trials,portion); 
end


%% Generate Event Structure

if findevnt_aud
    
    % Inputs for this function only
    soundpath = ['/Users/LauraLong/Documents/Lab/ECoGData/Task_Sounds/' task];
    audchannel = 1;
    confidencethreshold = .6;
    applyhilbert = [1:155];
    checkmanually = [];
    troubleshoot = [];
    stepsize = [];
    
    % Generate evnt structure!
    NeuralFindEvent_Audio_NAPLab(dpath,soundpath,subject,blocks,audchannel,confidencethreshold,applyhilbert,checkmanually,troubleshoot,stepsize);
end

if findevnt_trig
    
    % Inputs for this function only 
    soundpath = ['/Users/LauraLong/Documents/Lab/ECoGData/Task_Sounds/' task]; % directory containing .wav files 
    pulsepath = '/Users/LauraLong/Documents/Lab/ECoGData/NYUtrig.mat';  % path to original trigger pulse pattern (necessary if trigger is not already calculated to 0s/1s)
    audchannel = 1;
    troubleshoot = [];
    
    % Generage evnt structure!
    evnt = NeuralFindEvent_Trigger_NAPLab(dpath, soundpath, pulsepath, subject, blocks, audchannel, troubleshoot);
end


%% Generate Out Structure

if genout
    
    % Inputs for this function only
    befaft = [2 2];
    audchannels = 1;
    specflag = 'Auditory';
    switch whichdata
        case 'htkraw'
            outdataf = 400;
        otherwise
            outdataf = 100;
    end
    normout = 0;
    
    % Generate out structure!
    NeuralGenOut_NAPLab(dpath, whichdata, subject, blocks, channelinds, audchannels, outdataf, befaft, specflag, recordingsystem, normout);
end

%% Notify User Script Has Finished
disp('ModifyHTK Complete');
