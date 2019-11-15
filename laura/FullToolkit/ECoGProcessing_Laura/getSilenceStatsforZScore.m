%% Params
subject = '050_NY638';
silenceblock = 1;
saveblocks = 1;
dataorhtk = 'htk';
datacond = 'cleaned_highgamma';
basepath = ['/Users/LauraLong/Documents/Lab/ECoG Data/' subject filesep];

calcstats = 0;
normandsaveouts = 1;

savepath = [basepath 'processed' filesep 'B' num2str(saveblocks(i)) filesep 'silencestats.mat'];

if calcstats
    %% Load data
    
    switch dataorhtk
        case 'data'
            dpath = [basepath 'original' filesep 'B' num2str(silenceblock)];
            load([dpath filesep 'data.mat'])
        case 'htk'
            dpath = [basepath 'processed' filesep 'B' num2str(silenceblock) ];
            channels = getchannelinds([dpath filesep datacond], 'htk');
            data = [];
            for k = 1:length(channels)
                [data(k,:),fs_data] = readhtk([dpath filesep datacond filesep 'Ch' num2str(channels(k)) '.htk']);
            end
            audchannels = getchannelinds([dpath filesep 'analog'], 'htk');
            aud = [];
            for k = 1:length(audchannels)
                [aud(k,:),fs_aud] = readhtk([dpath filesep 'analog' filesep 'a' num2str(audchannels(k)) '.htk']);
            end
            
    end
    
    
    %% Plot all audio channels
    
    figure;
    for i = 1:size(aud,1)
        subplot(size(aud,1),1,i)
        plot(aud(i,:)); title(num2str(i));
    end
    suptitle('Audio Channels');
    
    
    %%
    bestaud = 1;
    figure; plot(aud(bestaud,:));
    
    %%
    
    recordblockend = 2.113e7; % end of landmark in aud channel (in samples)
    audioblockend = 8*60+20; % end of landmark in phone audio recording (in seconds)
    windowedges = [4*60+40 4*60+46; 8*60+10 8*60+17]; % the start and stop times of desired silence from phone audio recording (in seconds)
    
    recordblockend = recordblockend/fs_aud; % converts to seconds
    windowedges = round(windowedges - audioblockend + recordblockend); % adjust window edges to be the time after beginning of aud channel
    
    windowedges = [4.6e7 4.67e7; 4.69e7 4.76e7 ]/fs_aud;
    
    silwindow_data = []; silwindow_aud = [];
    for i = 1:size(windowedges,1)
        silwindow_data = round([silwindow_data windowedges(i,1)*fs_data:windowedges(i,2)*fs_data]);
        silwindow_aud = [silwindow_aud windowedges(i,1)*fs_aud:windowedges(i,2)*fs_aud];
    end
    
    silence_aud = aud(bestaud,silwindow_aud);
    figure; plot(silence_aud);
    
    silence_data = data(:,silwindow_data);
    figure; plot(silence_data');
    
    
    silmean = mean(silence_data,2);
    silstd = std(silence_data,[],2);
    
    figure;
    subplot(211)
    plot(silmean); title('Mean')
    subplot(212)
    plot(silstd); title('STD')
    
    %% save
    clear params;
    params.subject = subject;
    params.datacond = datacond;
    params.silenceblock = 1;
    params.silenceinds = windowedges;
    params.saveblocks = saveblocks;
    
    for i = 1:length(saveblocks)
        save(savepath,'silmean','silstd','params');
    end
    
end

%% if desired, run the zscore on an out structure

if normandsaveouts
    
    outfile = ['./outs/' subject '*.mat'];
    outfiles = dir(outfile);
    
    for i = 1:length(outfiles)
        thisoutfile = ['./outs' filesep outfiles(i).name];
        thisout = load(thisoutfile);
        out = zscore_out(thisout.out,savepath);
        save(thisoutfile,'out');
        disp(['saved ' thisoutfile])
    end
    
end


