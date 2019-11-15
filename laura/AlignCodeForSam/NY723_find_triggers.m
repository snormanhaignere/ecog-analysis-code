%% This script will process NYU data from .edf through evnt.mat
% Separate scripts written for each NYU patient; should have all data
% needed to replicate (timestamps, trigger channel, etc)

% Subject and folder info
subject = '086_NY723';
pth = ['/Users/LauraLong/Documents/Lab/ECoG Data/' subject filesep]; % SAM: Enter path to the location you have set up the proper file structure.

% Set up block info
block = 3; % SAM: Can change this if you have multiple blocks. 
% SAM: the below section would usually be what I use to select
% different edfs for different blocks; in your case you should only have
% one so can skip it. %
% switch block 
%     case 1
%         whichclin = 1; 
%     case 2
%         whichclin = 1;
%     case 3
%         whichclin = 2;
% end

%% Load data and plot to ID trigger channel
name = dir([pth 'original/all/*.EDF']);  % SAM: Change this to find the correct name of the edf file. 

% Load data
reload = 1;
if ~exist('hdr','var') || reload
    clear data_1 data_2 data_r data_s hdr*
    name = name.name; %     name = name(whichclin).name;
    [hdr,data_1] = edfread([pth 'original/all/' name]); % SAM: Change this to be path to edf file. Can take a long time to read. 
end


% SAM: Below section is a test to plot various channels so I can identify which one had the trigger. 
% In this case, this info is unlikely to change because I've already done this work for B1 and B2 of this patient.  
% However, if the triggers don't look correct, you can uncomment and use
% this section to find them. %
% Plot to check which is the trigger 
% disp(hdr.label');
% chantocheck = input('Which channels do you want to check? ');
% figure; 
% for i = 1:length(chantocheck)
%     subplot(length(chantocheck),1,i)
%     plot(data_1(chantocheck(i),:));
%     title(hdr.label(chantocheck(i)));
% end


%% Plot trigger to ID split times and plot blocks

% Plot full trigger to get split times
trigch = 165; % SAM: Here is where the identified trigger index goes. Unlikely to change in this case. 
figure; plot(data_1(trigch,:));

% Select split times and plot to check
switch block
    case 1
        blockrng = 5.242e4:8.266e5;
        task = 'Scrambling';
    case 2
        blockrng = 8.266e5:1.84e6;
        task = 'NYULocalizer';
    case 3
        blockrng = []; % SAM: Check out the trigger plot and identify the sample ranges that you want SAMPLE RANGE OF BLOCK HERE
        task = ''; % SAM: Name of folder where you will put task .wav files (names must match behavioral .mat file). 
end
figure; plot(data_1(trigch,blockrng));

% Split the data
alldata = data_1(:,blockrng);

%% Save data

% Identify channel ranges for each type of data
datachan = 1:164; % SAM: Enter range of data channels. Unlikely to change. 
trigch = 165; % SAM: Enter trigger channel index. Unlikely to change.
audch = 165; % SAM: Enter trigger channel index. Unlikely to change. (I double-save trigger as trig and aud when no aud is available)

% Split up the 
data = alldata(datachan,:); fs_data = hdr.samples(1)/hdr.duration; 
trig = alldata(trigch,:); fs_trig = fs_data;
aud = alldata(audch,:); fs_aud = fs_data;
chnames = hdr.label(datachan);
trignames = hdr.label(trigch);
orighdrs = hdr;

% Save data 
savedir = [pth 'original' filesep 'B' num2str(block)];
if ~exist(savedir,'dir'), mkdir(savedir),end
save([savedir '/data.mat' ],'data','trig','chnames','trignames','orighdrs','fs_trig','fs_data','aud','fs_aud','-v7.3');
save([savedir '/allchnames_B' num2str(block) '.mat'],'chnames');


%% Convert to HTK to get trigger files ready for processing

ConvertToHTK_Laura(pth,'edf',block,[],[]);

%% Get and write evnt-friendly triggers

trigpath = [pth 'processed/B' num2str(block) '/trigger/'];

[t1,fs1] = readhtk([trigpath 't1.htk']); 
a = load('/Users/LauraLong/Documents/Lab/ECoG Data/NYUtrig.mat');  % SAM: Change this to the location of the NYUtrig file I send you. 
t2 = a.trigform; fs2 = a.fs;
[trigger,fs_trig] = matchTrigger(t1,fs1,t2,fs2);

save([trigpath 'trigger.mat'],'trigger','fs_trig');
writehtk([trigpath 'trigger.htk'],trigger,fs_trig);

%% Now get stimulus and evnt info

stimpath = '/Users/LauraLong/Documents/Lab/ECoG Data/Task_Sounds'; % SAM: Change this to the path where you have the stimulus folder. 
soundpath = [stimpath filesep task]; % SAM: Change this to the path where you have the .wav files of the stimuli. 

makeStimFolder_Laura(pth,block,stimpath,task);

trigchannel = 1;
troubleshoot = 0; 
evnt = NeuralFindEvent_Trigger(pth, soundpath, subject, block, trigchannel);


%%

rates = [sqrt(1/3), sqrt(3)];



