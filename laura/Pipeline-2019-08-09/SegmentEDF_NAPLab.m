%% This script will segment .edfs into separate blocks
% This is most commonly necessary for NYU data, and must be done before the steps in ModifyHTK.
% The edf you desire to split should be saved in pth/original/all
% Blocks of code that user must change are denoted with 'USER:'

% USER: Subject and folder info
subject = '086_NY723'; % subject name 
pth = ['/Users/LauraLong/Documents/Lab/ECoGData/' subject filesep]; % Path to the location you have set up the proper file structure
block = 3; % which block you are currently processing; run this script once per block

%% Load data
name = dir([pth 'original/all/*.edf']); 

% Load data
reload = 0; 
if ~exist('hdr','var') || reload
    clear data_1 data_2 data_r data_s hdr*
    name = name.name; %     name = name(whichclin).name;
    [hdr,data_1] = edfread([pth 'original/all/' name]); % Note: for large files, this line can take a very long time!
end

%% Display channel names and plot data to ID aud, trig, and neural channels

% Display labels
alllabels = hdr.label; numlabels = length(alllabels)
for i = 1:numlabels
    disp([num2str(i) '     ' alllabels{i}]);
end

% Plot subsets of channels to ID channels
chantocheck = input('Which channels do you want to check? '); % any channels you want to see plotted
if ~isempty(chantocheck)
    figure;
    for i = 1:length(chantocheck)
        subplot(length(chantocheck),1,i)
        plot(data_1(chantocheck(i),:));
        title(hdr.label(chantocheck(i)));
    end
end

% USER: Enter channel ranges for each type of data
datachan = 1:164; % Enter range of data channels
trigch = 165; % Enter trigger channel index
audch = 165; % Enter audio channel index (if no audio is available, I double-save trigger as audio)


%% Identify block split times

% Plot trigger channel so user can find samples
figure; plot(data_1(trigch,:));

% USER: Check out the trigger plot and identify the sample ranges that you want 
% Note: it's recommended to leave a bit of extra time at the front/end 
% (the out structure may include a couple seconds before and after each stim)
switch block
    case 1
        blockrng = 5.242e4:8.266e5; % SAMPLE RANGE OF BLOCK 1 HERE
    case 2
        blockrng = 8.266e5:1.84e6; % SAMPLE RANGE OF BLOCK 2 HERE
    case 3 % ... and so forth; add more cases for additional blocks
        blockrng = []; 
end

% Plot split to confirm
figure; plot(data_1(trigch,blockrng));


%% Split and save data

% Grab just this block range
alldata = data_1(:,blockrng);

% Now split into data, trig, aud
data = alldata(datachan,:); fs_data = hdr.samples(datachan(1))/hdr.duration; chnames = hdr.label(datachan);
trig = alldata(trigch,:); fs_trig = hdr.samples(trigchan(1))/hdr.duration; trignames = hdr.label(trigch);
aud = alldata(audch,:); fs_aud = hdr.samples(audchan(1))/hdr.duration; audnames = hdr.label(audch);
orighdrs = hdr;

% Save data 
savedir = [pth 'original' filesep 'B' num2str(block)];
if ~exist(savedir,'dir'), mkdir(savedir),end
save([savedir '/data.mat' ],'data','fs_data','chnames','trig','fs_trig','trignames','aud','fs_aud','audnames','orighdrs','-v7.3');
save([savedir '/allchnames_B' num2str(block) '.mat'],'chnames');
