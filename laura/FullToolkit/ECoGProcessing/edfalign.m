%% Analyze Lucia's Data
subject = '050_NY638';
block = 1;
basepath = ['/Users/LauraLong/Documents/Lab/ECoG Data/' subject filesep];
edfpath = [basepath 'original' filesep 'B' num2str(block) filesep];
loaddata = 0;
getinput = 0;

if loaddata
    fnames = dir([edfpath '*.edf']);
    hdr = cell(1, length(fnames)); record = hdr;
    if getinput, trigch = hdr; audch = hdr; datach = hdr; end
    ecog = hdr; rawtrig = hdr; aud = hdr;
    for i = 1:length(fnames)
        [hdr{i}, record{i}] = edfread([edfpath fnames(i).name]);
        
        if getinput
            
            if isfield(hdr{i},'label')
                for j = 1:length(hdr{i}.label)
                    disp([num2str(j) hdr{i}.label(j)]);
                end
            end

            datach{i} = input('Which channels are data? ');
            try allfs_data = hdr{i}.samples(datach{i});
                allfs_data = unique(allfs_data);
                if length(allfs_data) == 1
                    fs_data = allfs_data;
                else
                    error('Multiple fs for data channels');
                end
            catch
                fs_data = input('What is fs_data? ');
            end
            trigch{i} = input('Which channels are trigger? ');
            try allfs_trig = hdr{i}.samples(trigch{i});
                allfs_trig = unique(allfs_trig);
                if length(allfs_trig) == 1
                    fs_trig = allfs_trig;
                else
                    error('Multiple fs for data channels');
                end
            catch
                fs_trig = input('What is fs_trig? ');
            end
            audch{i} = input('Which channels are audio? ');
            try allfs_aud = hdr{i}.samples(audch{i});
                allfs_aud = unique(allfs_aud);
                if length(allfs_aud) == 1
                    fs_aud = allfs_aud;
                else
                    error('Multiple fs for data channels');
                end
            catch
                fs_aud = input('What is fs_aud? ');
            end
        end
        
        allchnames{i} = hdr{i}.label(datach{i});
        ecog{i} = record{i}(datach{i},:);
        rawtrig{i} = record{i}(trigch{i},:);
        audio{i} = record{i}(audch{i},:);
        
    end
end

%% Process Triggers

trig = cell(1,length(rawtrig));
triglocs = trig; if getinput, goodtrigch = trig; end
for i = 1:length(rawtrig)
        
    if getinput
        figure;
        for j = 1:size(rawtrig{i},1)
            subplot(size(rawtrig{i},1),1,j);
            plot(rawtrig{i}(j,:));
            title(num2str(j));
        end
        goodtrigch{i} = input('Which trigger channels are good? ');
        rawtrig{i} = rawtrig{i}(goodtrigch{i},:);
    end
    
    rawtrig{i} = sum(rawtrig{i},1); % sum all channels together
    trig{i} = rawtrig{i} > max(rawtrig{i})*.5; % find samples where trigger is at least 80%
    trig{i} = [false, diff(trig{i}) > 0];
    triglocs{i} = find(trig{i});

    

    %% ADJUST TRIGGERS IF NECESSARY
    % Plot original
    figure; subplot(211); plot(trig{i}); title(['original trigger: ' num2str(sum(trig{i}))])
    
    % Deletions
    switch i 
        case 1
            trig{i}(triglocs{i}(1:2)) = 0;
        case 2
    end

    % Thresholding in time
    trig{i}(triglocs{i}(triglocs{i}>14e5)) = 0;
    
    % Replot
    subplot(212); plot(trig{i}); title(['revised trigger: ' num2str(sum(trig{i}))])
    suptitle(['record ' num2str(i)]);
    
    % Refind locations
    triglocs{i} = find(trig{i});
    
    
end


%% Check triggers with audio

% Find the audio files
for i = 1:length(audio)
    
    if getinput && size(audio{i},1) > 1
        figure;
        for j = 1:size(audio{i},1)
            subplot(size(audio{i},1),1,j);
            plot(audio{i}(j,:)); title(j)
        end
        audchannel = input('Which audio channel is best? Input ''0'' for none: '); % INPUT GOOD TRIGGER CHANNELS HERE; alternatively audchannel = [];
    else
        audchannel = 1;
    end
    audio{i} = audio{i}(audchannel,:);
    
    % Plot audio and triggers together to compare and look for mistakes
    figure;
    if fs_aud == fs_trig
        plot(trig{i}*1e4,'k','LineWidth',2); hold on; 
        title(['num triggers: ' num2str(sum(trig{i}))]);
        plot(audio{i}); axis tight
    else
        figure;
        subplot(211);
        plot(trig{i}); axis tight
        title(['num triggers: ' num2str(sum(trig{i}))]);
        subplot(212);
        plot(audio{i}); axis tight
        title('audio');
    end

end



%% Align edf 1 and 2

% Figure out min and max of triggers
trigbegin = min([triglocs{:}]) - 1;
trigend = min( length(trig{1}) - triglocs{1}(end) , length(trig{2}) - triglocs{2}(end));

align = cell(1,length(trig)); trigalign = align;
for i = 1:length(trig)
    align{i} = [triglocs{i}(1) - trigbegin  : triglocs{i}(end) + trigend];
    trigalign{i} = trig{i}(align{i});
    % Adjustments if required
    switch i
        case 1
            
        case 2
            trigalign{i} = [0 0 0 trigalign{i}];
    end
end

minlength = min(length(align{1}),length(align{2}));
for i = 1:length(trigalign)
    trigalign{i} = trigalign{i}(1:minlength);
    align{i} = align{i}(1:minlength);
end

figure; subplot(211), plot(trigalign{1}); hold on; plot(trigalign{end}); hold off; title('Triggers Overlaid');
subplot(212), plot(find(trigalign{1}) - find(trigalign{2})); title('Trigger Indices Subtracted (1-2)'); addonsetline(0,'horizontal');

% Determine trigger to save
if getinput
    besttrigch = input('Which trigger is better? '); 
    trigger = trigalign{besttrigch};
else
    trigger = trigalign{i};
end


%% Now clip and save the data according to the new alignment indices

aud = []; data = []; chnames = [];
for i = 1:length(trigalign)
    data = cat(1,data,ecog{i}(:,align{i}));
    aud = cat(1,aud,audio{i}(:,align{i}));
    chnames = cat(2,chnames,allchnames{i});
end



%% Save files

save([edfpath filesep 'data.mat'],'data','aud','trigger','fs_data','fs_aud','fs_trig');
save([edfpath filesep 'allchnames_B' num2str(block) '.mat'],'chnames');


triggerpath = [basepath filesep 'processed/B' num2str(block) filesep 'trigger'];
if ~exist('triggerpath','dir'), mkdir(triggerpath); end
save([triggerpath filesep 'trigger.mat'],'trigger','fs_trig');
writehtk([triggerpath filesep 'trigger.htk'],trigger,fs_trig);

