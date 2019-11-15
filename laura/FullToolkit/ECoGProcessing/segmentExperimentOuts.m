% Segment experiment into separate outs. Assumes evnt has already been written.

getinput = 1;

%% Setup
subject = '050_NY638';
task = 'Segmentation';
dpath = ['/Users/LauraLong/Documents/Lab/ECoG Data/' subject];
recordingsystem = 'edf';

% Enter which data you want to start with.
% If empty, defaults to htkraw
whichdata = 'cleaned_highgamma';

% Enter additional parameters.
% If empty, defaults to all blocks, all channels, and 1000 Hz.
block =  1;
channels = [];
    

%% Correct to 
switch lower(task)
    case 'segmentation'
%         [ phraseevnt, testoldevnt ] = segmentationEvntTesttoPhrase(subject,blocks,fixstimorevnt,keeporig )
end


%%


load([dpath filesep 'processed' filesep 'B' num2str(block) filesep 'evnt_' subject, '_B' num2str(block) '.mat']);
if getinput
    for i = 1:length(evnt)
        disp([num2str(i) '        ' evnt(i).name]);
    end
    breakpoints = input('Where are the breakpoints? ');
    evntnames = input('What are the evntnames? ');
else
    breakpoints = [25 62 86];
    evntnames = {'pre' 'train' 'post' 'control'};
end


%% Find Windows

windows = cell(1,length(breakpoints)+1);
start = 1;
for i = 1:length(windows)
    if i == length(windows)
        windows{i} = start:length(evnt);
    else
        windows{i} = start:breakpoints(i)-1;
        start = breakpoints(i);
    end
    disp(['New Evnt: ' evntnames{i}]);
    disp({evnt(windows{i}).name}');
end

%% Generate Out Structure

befaft = []; % before/after window for each stimulus; defaults to [.5 .5]
audchannels = 1; % which aud channels to include from the recording
specflag = 'Auditory'; % type of spectrogram
if isempty(whichdata) % fs of data defaults to 400 for raw, 100 otherwise
    outdataf = 400;
    whichdata = 'htkraw';
else
    outdataf = 100;
end


% Loop over evnts and save separate outs
for i = 1:length(windows)
    out = NeuralGenOut_Laura(dpath, whichdata, subject, block, channels, audchannels, outdataf, befaft, specflag, recordingsystem, windows{i});
    save([dpath filesep 'processed' filesep 'B' num2str(block) filesep 'out_' subject, '_B' num2str(block) '_' whichdata '_' evntnames{i} '.mat'],'out','-v7.3');
    
end

