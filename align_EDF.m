function align_EDF(exp, subjid, block, varargin)

% exp = 'speech-TCI';
% subjid = '086_NY723';
% block = 3;

global root_directory;
project_directory = [root_directory '/' exp];
ecog_raw_directory = [project_directory  '/data/ECoG-Raw/' subjid '/'];

%% Optional parameters

clear I;
I.edf_sr = 1024; % sampling rate of EDF data
I.trigchan = NaN; 
I.blockrng = []; % range of samples for this block
I.ecogchan = []; % channels with ECoG data
I.audchan = NaN;
I.trigMATfile = [project_directory '/analysis/trig.mat'];
I.overwrite = false;
I = parse_optInputs_keyvalue(varargin, I);

%% Create data file

edf_directory = [project_directory '/data/ECoG-EDF/' subjid '/']; % SAM: Enter path to the location you have set up the proper file structure.
save_directory = [edf_directory 'original/B' num2str(block) '/'];
data_MAT_file = [save_directory 'data.mat'];
channel_MAT_file = [save_directory 'allchnames_B' num2str(block) '.mat'];

if ~exist(data_MAT_file, 'file') || ~exist(channel_MAT_file, 'file') || I.overwrite
    
    %% Load EDF file
    
    % path to EDF directory
    edf_file = [edf_directory 'original/B' num2str(block) '/data_' num2str(I.edf_sr) '.EDF'];
    assert(logical(exist(edf_file,'file')));
    assert(eq_tol(hdr.samples(1)/hdr.duration, I.edf_sr));
    [hdr,data_1] = edfread(edf_file);
    
    %% Find trigger channel
    if isnan(I.trigchan)
        figh = figure;
        while true
            for j = 1:length(hdr.label)
                fprintf('%d: %s\n',j,hdr.label{j});
            end
            chantocheck = input('Channels to check (-1 when done): ');
            clf(figh);
            for i = 1:length(chantocheck)
                subplot(length(chantocheck),1,i)
                plot(data_1(chantocheck(i),:));
                title(sprintf('chan: %d, %s', chantocheck(i), hdr.label{chantocheck(i)}));
            end
            I.trigchan = input('Trigger channel (NaN if not found): ');
            if ~isnan(I.trigchan)
                break;
            end
        end
    end
    
    %% Choose range of samples
    
    if isempty(I.blockrng)
        figure;
        plot(data_1(I.trigchan,:));
        I.blockrng = input('Block range in samples: ');
    end
    
    % Select samples from this block range
    alldata = data_1(:,I.blockrng);
    
    %% Save data
    
    % Split up the
    data = alldata(I.ecogchan,:);
    fs_data = hdr.samples(1)/hdr.duration;
    fs_trig = fs_data;
    fs_aud = fs_data;
    if ~isempty(I.trigchan)
        trig = alldata(I.trigchan,:);
        trignames = hdr.label(I.trigchan);
    else
        trig = [];
        trignames = {};
    end
    if ~isempty(I.audchan)
        aud = alldata(I.audchan,:);
    else
        aud = [];
    end
    chnames = hdr.label(I.ecogchan);
    orighdrs = hdr;
    
    % Save data
    if ~exist(save_directory,'dir'); mkdir(save_directory); end
    save(data_MAT_file,'data','trig','chnames','trignames','orighdrs','fs_trig','fs_data','aud','fs_aud','-v7.3');
    save(channel_MAT_file,'chnames');
    
else
    
    load(data_MAT_file);
    load(channel_MAT_file);
end

%% Convert to HTK to get trigger files ready for processing

ConvertToHTK_Laura(edf_directory,'edf',block,[],[]);

%% Get and write evnt-friendly triggers

trigpath = [edf_directory 'processed/B' num2str(block) '/trigger/'];

[t1,fs1] = readhtk([trigpath 't1.htk']); 
assert(logical(exist(I.trigMATfile,'file')));
a = load(I.trigMATfile);  % SAM: Change this to the location of the NYUtrig file I send you. 
t2 = a.trigform; fs2 = a.fs;
[trigger,fs_trig] = matchTrigger(t1,fs1,t2,fs2);

save([trigpath 'trigger.mat'],'trigger','fs_trig');
writehtk([trigpath 'trigger.htk'],trigger,fs_trig);

%% Now get stimulus and evnt info

stimpath = [project_directory '/stimuli/stimorder/' subjid];
soundpath = [project_directory '/stimuli/final-stims']; % SAM: Change this to the path where you have the .wav files of the stimuli. 
task = '';  
makeStimFolder_Laura(edf_directory,block,stimpath,task);

trigchannel = 1;
troubleshoot = 0; 
evnt = NeuralFindEvent_Trigger(edf_directory, soundpath, subjid, block, trigchannel, troubleshoot);

%% create a link to the appropriate directory to make compatible with the rest of Sam's pipeline

processed_directory = [edf_directory 'processed/B' num2str(block)];
link_directory = mkpdir([ecog_raw_directory 'B' num2str(block)]);
if ~exist(link_directory, 'dir')
    unix(['ln -s ' processed_directory ' ' link_directory]);
end
copyfile(channel_MAT_file, [ecog_raw_directory 'allchnames_B1.mat']);

