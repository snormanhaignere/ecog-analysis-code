function ConvertBlackRockToHTK 

%Set variables (make these inputs later, but for debugging I'm declaring them)
audiomatfile = 'C:\Users\Laura\Documents\MATLAB\NimaLab\20150205-114247\20150205-114247-001 - ns5.mat';
neuralmatfile = 'C:\Users\Laura\Documents\MATLAB\NimaLab\20150205-114247\20150205-114247-001 - ns4.mat';
NumOfBlocks = 4;
destpath = 'C:\Users\Laura\Documents\MATLAB\NimaLab\20150205-114247';


audiomat = load(audiomatfile);
audiodata = double(audiomat.(['C' num2str(audiomat.Header.ChannelID)]));
soundf = audiomat.Header.Fs;

neuralmat = load(neuralmatfile);
numchannels = neuralmat.Header.ChannelCount;
neuraldata = zeros(numchannels, length(double(neuralmat.(['C' num2str(neuralmat.Header.ChannelID(1))]))) );
for cnt1 = 1:numchannels
    neuraldata(cnt1,:) = double(neuralmat.(['C' num2str(neuralmat.Header.ChannelID(cnt1))]));
    display('loaded channel');
end
dataf = neuralmat.Header.Fs;

% convert interop continous recording to standard blocks

% display the audio and ask the person to specify blocks
plot(audiodata(1:100:end));
tmp = ginput(NumOfBlocks+1);
% convert to sample time:
tmp = 100*tmp/soundf;
tmp = tmp(:,1);

for cnt1 = 1:NumOfBlocks-1
    % cut the sound and neural data for each block
    aud_tmp  = audiodata(tmp(cnt1)*soundf:tmp(cnt1+1)*soundf);
    data_tmp = neuraldata(tmp(cnt1)*dataf:tmp(cnt1+1)*dataf);
    % where to write this?
    datapath = [destpath filesep 'B' num2str(cnt1) filesep 'htkraw' filesep];
    if ~exist(datapath,'dir'), mkdir(datapath);end
    stimpath = [destpath filesep 'B' num2str(cnt1) filesep 'analog' filesep];
    if ~exist(stimpath,'dir'), mkdir(stimpath);end
    % now write the sound of this block to the "aoppropriate directory
    writehtk([stimpath filesep 'a1.htk'],aud_tmp,soundf);
    for cnt2 = 1:size(neuraldata,1)
        htkname = [datapath 'Ch' num2str(cnt2) '.htk'];
        writehtk(htkname,datapath(cnt2,:),dataf);
    end
end
    
