function ConvertBlackRockToHTK(audiomatfile,neuralmatfile,NumOfBlocks,destpath)
% convert interop continous recording to standard blocks

% load the audio
audiomat = load(audiomatfile);
audiodata = double(audiomat.(['C' num2str(audiomat.Header.ChannelID)]));
soundf = audiomat.Header.Fs;

% display the sound and ask the subject to segment
plot(audiodata(1:100:end));
tmp = ginput(NumOfBlocks+1);
% convert to sample time:
tmp = 100*tmp/soundf;
tmp = tmp(:,1);

%load the neural data
neuralmat = load(neuralmatfile);
dataf = neuralmat.Header.Fs;
numchannels = neuralmat.Header.ChannelCount;
neuraldata = zeros(1,length(double(neuralmat.(['C' num2str(neuralmat.Header.ChannelID(1))]))) );
for cnt1 = 1:numchannels % channel number
    neuraldata = double(neuralmat.(['C' num2str(neuralmat.Header.ChannelID(cnt1))]));
    display('loaded channel');
    for cnt2 = 1:NumOfBlocks % block number
        % cut the sound and neural data for each block
        aud_tmp  = audiodata(tmp(cnt2)*soundf:tmp(cnt2+1)*soundf);
        data_tmp = neuraldata(tmp(cnt2)*dataf:tmp(cnt2+1)*dataf);
        % where to write this?
        datapath = [destpath filesep 'B' num2str(cnt2) filesep 'htkraw' filesep];
        if ~exist(datapath,'dir'), mkdir(datapath);end
        stimpath = [destpath filesep 'B' num2str(cnt2) filesep 'analog' filesep];
        if ~exist(stimpath,'dir'), mkdir(stimpath);end
        % now write the sound of this block to the "aoppropriate directory
        writehtk([stimpath filesep 'a1.htk'],aud_tmp,soundf);
        htkname = [datapath 'Ch' num2str(cnt1) '.htk'];
        writehtk(htkname,data_tmp,dataf);
    end
end

