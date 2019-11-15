%% User Inputs: Set parameters, choose what actions to take

% Enter subject ID.
subject = '085_CU102';
task = 'NYULocalizer';

% Enter dpath. Should be the patient folder, containing 'original' directory.
% If empty, defaults to current directory.
dpath = ['/Users/LauraLong/Documents/Lab/ECoG Data/' subject];

% Enter recording system. Necessary for ConvertToHTK, otherwise optional.
% Options: 'tdt' 'blackrock' 'xltek 'natus' 'eeg' 'newtdt'
recordingsystem = 'blackrock';

% Enter which data you want to start with.
% If empty, defaults to htkraw.
whichdata = '';

% Enter additional parameters.
% If empty, defaults to all blocks, all channels, and 1000 Hz.
blocks = [2]; 
channelinds = [];
dataf = [];

% Possible data processing steps. Set to 1 to take each action.
% Will execute in order listed.
converttohtk = 0;
writechnames = 0;
findcleanvars = 0; cleanblock = 1; blocks(1); % select the desired block to use to select cleaning vars
cleandata = 0; 
highgamma = 0;
freqband = 0;
makestim = 1;
findevnt = 1;
findevnt_trig = 0;
genout = 0; normout = 1;

plotresults = 1; % say yes if you want a ton of figures

% Extra params are necessary for CMR, stim, evnt, and out. 
% If applicable, update them in sections below!srgsr


%%  Convert to HTK
if converttohtk
    ConvertToHTK_Laura(dpath,recordingsystem,blocks,channelinds,dataf);
end

if writechnames
    trodenames = {};
    trodelens = [];
    WriteChNames_Laura(dpath,blocks,trodenames,trodelens);
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
    stimpath = '/Users/LauraLong/Documents/Lab/ECoG Data/Task_Sounds';
    trials = [1];
    portion = {[1:4 6:19 20 145:155]};
%     portion = {5};
%     portion = {[1:36 38:144]};
%     portion = {[1:48 50:67 69]}; % for localizer
%     portion = {[]};
     
    % Generate stim folders!
    makeStimFolder_Laura(dpath,blocks,stimpath,task,trials,portion);
end


%% Generate Event Structure

if findevnt
    
    % Inputs for this function only
    soundpath = ['/Users/LauraLong/Documents/Lab/ECoG Data/Task_Sounds/' task];
    audchannel = 1;
    troubleshoot = 0;
    needhilb = [1:155];
    confidencethreshold = .1;
    checkmanually = []; 
    
    
    % Generate evnt structure!
    NeuralFindEvent_Hilbert(dpath,soundpath,subject,blocks,audchannel,troubleshoot,needhilb,confidencethreshold,checkmanually);
end

if findevnt_trig
    audchannel = 1;
    soundpath = ['/Users/LauraLong/Documents/Lab/ECoG Data/Task_Sounds/' task];
    troubleshoot = 0;
    evnt = NeuralFindEvent_Trigger(dpath, soundpath, subject, blocks, audchannel, troubleshoot);
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
    
    % Generate out structure!
    NeuralGenOut_Laura(dpath, whichdata, subject, blocks, channelinds, audchannels, outdataf, befaft, specflag, recordingsystem, normout);
end

%% Notify User Script Has Finished
disp('ModifyHTK Complete');
