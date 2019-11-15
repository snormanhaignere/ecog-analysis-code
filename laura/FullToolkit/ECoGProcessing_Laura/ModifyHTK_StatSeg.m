%% Statistics and Segmentation Specifically
%% User Inputs: Set parameters, choose what actions to take

% Enter subject ID.
subject = '055_CUBF35'
task = 'Localizer';

% Enter relevant paths. Dpath should be the patient folder, containing 'original' directory.
% If empty, defaults to current directory.
dpath = ['/Users/LauraLong/Documents/Lab/ECoG Data/' subject];
stimpath = '/Users/LauraLong/Documents/Lab/ECoG Data/Task_Sounds';
soundpath = ['/Users/LauraLong/Documents/Lab/ECoG Data/Task_Sounds/' task];

% Enter recording system. Necessary for ConvertToHTK, otherwise optional.
% Options: 'tdt' 'blackrock' 'xltek 'natus' 'eeg' 'newtdt'
recordingsystem = 'edf';

% Enter which data you want to start with.
% If empty, defaults to htkraw.
whichdata = '';

% Enter additional parameters.
% If empty, defaults to all blocks, all channels, and 1000 Hz.
blocks = [2];
trials = [1];
channelinds = [];
dataf = [];

% Possible data processing steps. Set to 1 to take each action.
% Will execute in order listed.
converttohtk = 0;
findcleanvars = 0; cleanblock = 1; % select the desired block to use to select cleaning vars
cleandata = 1; 
highgamma = 1;
freqband = 0;
makestim = 0;
gettrig = 0;
findevnt = 0;
findevnt_trig = 0;
getsilstat = 0;
genout = 1; normout = 1;
savefft = 0;
saveunits = 0;

plotresults = 1; % say yes if you want a ton of figures

% Extra params are necessary for CMR, stim, evnt, and out. 
% If applicable, update them in sections below!


%%  Convert to HTK (chnames happens here)
if converttohtk
    ConvertToHTK_Laura(dpath,recordingsystem,blocks,channelinds,dataf);
end

%% Clean data


if findcleanvars
    cleanvars = findCleanDataVariables(dpath,cleanblock,task);
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
    
%     filters.highpassfreq = .1;
    
    % Clean Data!
    for i = 1:length(blocks)
        HTKtoCleanData(dpath,whichdata,blocks(i),channelinds,dataf,filters.highpassfreq,reftype,amplifiers,channels.bad,epochs,channels.iffy,normdata,plotresults);
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
    
    HTKtoHighGamma_Laura(dpath,whichdata,blocks,channelinds,dataf);
    
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
    
    HTKtoFreqBand(dpath,whichdata,blocks,channelinds,dataf,freqRange,getenvelope,datatag);
    
    if isempty(whichdata)
        whichdata = datatag;
    else
        whichdata = [whichdata '_' datatag];
    end
end

%% Make Stimulus Folders

if makestim
    
    % Inputs for this function only
    portion = {};
     
    % Generate stim folders!
    makeStimFolder_Laura(dpath,blocks,stimpath,task,trials,portion);
end


%% Get Triggers If Desired (For Statistical- gets Syllable/Phrase Timing and Trigger Vector)

if gettrig && strcmpi(task,'statistical')
    
    % Input for this function only
    trigchannel = 1;
    
    % Get NYU triggers!
    getNYUTriggers(dpath,blocks,trigchannel);
    
end

%% Generate Event Structure

if findevnt
    
    % Inputs for this function only
    audchannel = 1;
    troubleshoot = 0;
    needhilb = [1:201];
    confidencethreshold = .25;
    
    % Generate evnt structure!
    NeuralFindEvent_Hilbert(dpath,soundpath,subject,blocks,audchannel,troubleshoot,needhilb,confidencethreshold);
end

if findevnt_trig
    audchannel = 1;
    troubleshoot = 0;
    evnt = NeuralFindEvent_Trigger(dpath, soundpath, subject, blocks, audchannel, troubleshoot);
end

%% Generate Out Structure

if genout
    
    % Inputs for this function only
    befaft = [];
    audchannels = 1;
    specflag = 'Auditory';
    switch whichdata
        case 'htkraw'
            outdataf = 400;
        otherwise
            outdataf = 100;
    end
    
    % Generate out structure!
    NeuralGenOut_Laura(dpath, whichdata, subject, blocks, channelinds, audchannels, outdataf, befaft, specflag, recordingsystem, normout);
    
    % Update data tag 
    if isempty(whichdata)
        whichdata = 'htkraw_norm';
    else
        whichdata = [whichdata '_norm'];
    end
end


%% Save Statistics

if savefft
    
    savepath = dpath; % if savepath is dpath, it will save to the block folder
    fftlength = 250;
    hamminglength = 50;
    usefftnumfreqs = 1;
    
    saveFFTData(dpath,savepath,whichdata,subject,blocks,trials,fftlength,hamminglength,usefftnumfreqs,plotresults);
end

if saveunits
    
    savepath = dpath; % if savepath is dpath, it will save to the block folder
    alignstartorstop = 'start';
    window = [];
    compressresp = 1;
    arrangebymanner = 0;
    outdataf = 100;

    fnames = {'aud' 'sound' }; %'presp' 'presp_norm' 'rstim_all'}; 'rstim_syl'; 'rstim_wrd';
    
    switch task
        case 'Statistical'
            saveStatUnitData(dpath,soundpath,savepath,whichdata,subject,blocks,trials,fnames,outdataf,window,alignstartorstop,compressresp,normout,plotresults,arrangebymanner);
        case 'Segmentation'
    end
    
end


%% Notify User Script Has Finished
disp('ModifyHTK Complete');
