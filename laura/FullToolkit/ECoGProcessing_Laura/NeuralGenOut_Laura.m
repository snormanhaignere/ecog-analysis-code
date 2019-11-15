function [out] = NeuralGenOut_Laura(datapath, cond, subject, blocks, channels, audchannels, dataf, befaft, specflag, datatype, normout, evntpart, artifact)
% [out] = NeuralGenOut_Laura(datapath, cond, subject, blocks, channels, dataf, befaft, specflag, datatype, artifact)
%
% Inputs:
%       datapath (optional): subject folder
%           datapath must contain 'Stimulus' folder containing StimOrder.mat
%           within its 'processed' subdirectory
%               eg: datapath / processed / B# / Stimulus
%           default: current directory
%       cond: which data to use in out
%           options: directory names within processed block folders
%           default: 'htkraw'
%       subject: subject ID number (will be used to name out)
%       blocks (optional): specify which blocks to process
%           format: array of integers OR cell array of block folder names
%               eg: [1:4] OR {'B1' 'B2' 'B3' 'B4'}
%           default: all blocks in datapath / original
%       channels (optional): specify which channels to process
%           default: to all channels in each block
%       dataf (optional): specify desired data fs; best if a multiple of 1000
%           default: 1000Hz
%       befaft (optional): amount of silence to insert before /after each evnt
%       specflag: what type of spectrogram
%           options: Auditory, AudNorm, mfcc
%           default: Auditory
%       datatype (optional): specify data format
%           options: 'tdt' 'blackrock' 'xltek' 'natus' 'eeg'
%           default: 'ECoG'
%       artifact (optional): artifacts that you want to add
%           must be included in artifact folder
%           default: empty
%
% Result:
%       For each block, reads evnt and aligns to ECoG data
%       Clips ECoG data
%       Generates spectrogram
%       Asks user whether to save out file
%
%       Resulting file:
%           Out: datapath / processed / B# / out_subject_cond_B#.mat
%
% This version adapted from Bahar's code by Laura 11/2015;
% - removed evnt input; instead, it takes a blocks input & loads
% the appropriate evnt acccordingly. Allows looping over several blocks at
% once, and is more consistent with evnt code inputs.
% - adjusted for new ECoG Data organization (processed folder)
% - added option to save out if desired, saves to each block folder
% - added lots of comments
% - revised blocks input to be valid as matrix or cell
% - changed variable names and order (to lowercase, switched input order to be
% consistent with Convert/Out functions)

% Last update 12/10/15 by Laura
% Moved spectrogram to after the dataload to fix dataf problem

%% Deal with Inputs

fprintf('\nNeuralGenOut\n');

% Check if datapath was entered; if not, default to current directory
if ~exist('datapath','var') || isempty(datapath)
    datapath = pwd;
end

% Check if cond was entered; if not, default to htkraw
if ~exist('cond','var') || isempty(cond)
    cond = 'htkraw';
end

% Check if blocks was entered; if not, default to
if ~exist('blocks','var') || isempty(blocks)
    blocks = [];
end
if iscell(blocks) % If blocks input is matrix, turn it into a cell array
    blocknames = blocks;
else
    for i = 1:length(blocks)
        blocknames{i} = ['B' num2str(blocks(i))];
    end
end

% Check if channels was entered; if not, default to .5s on both sides
if ~exist('channels','var') || isempty(channels)
    channels = getchannelinds([datapath filesep 'processed' filesep blocknames{1} filesep cond],'htk'); % find channel inds
end

% Check if dataf was entered; if not, set resample flag to 0
if ~exist('dataf','var') || isempty(dataf)
    no_resample = 1;
else
    no_resample = 0;
end

% Check if befaft was entered; if not, default to .5s on both sides
if ~exist('befaft','var') || isempty(befaft)
    befaft=[0.5,0.5];
end

% Check if datatype was entered; if not, default to ECoG
if ~exist('datatype','var') || isempty(datatype)
    datatype = 'ECoG';
end

% Check if datatype was entered; if not, default to ECoG
if ~exist('specflag','var') || isempty(specflag)
    specflag = 'Auditory';
end

% Check if artifact was entered; if not, default to empty
if ~exist('normout','var') || isempty(normout)
    normout=0;
end


% Check if artifact was entered; if not, default to empty
if ~exist('artifact','var') || isempty(artifact)
    artifact=[];
end


% Check if evntpart was entered; if not, default to full evnt
if ~exist('evntpart','var') || isempty(evntpart)
    fullevnt = 1;
else
    fullevnt = 0;
end

%% Loop over blocks

for i = 1:length(blocknames)
    
    disp(['Processing ' blocknames{i} '...']);
    
    % Reset evnt, out, data; load the block's evnt
    evnt = []; out = []; record = [];
    load([datapath filesep 'processed' filesep blocknames{i} filesep 'evnt_' subject, '_' blocknames{i}, '.mat']);
    if ~fullevnt
        evnt = evnt(evntpart);
    end
    
    % Load evntnames, make out.name appropriate length
    evntnames = {evnt.name}';
    [~,unique_mat,unique_index]=unique(evntnames); % find unique names
    out = struct('name',cell(1,length(unique_mat)));
    
    prevblock='  '; % previous block; keeps track of whether the block has changed and new data should be loaded; starts blank
    loadload; close; % loads audio colormap, filterbank, params, gets rid of resulting figure
    
    % Load chnames
    if exist([datapath filesep 'original' filesep blocknames{i} filesep 'allchnames_' blocknames{i} '.mat'],'file')
        allchnames = load([datapath filesep 'original' filesep blocknames{i} filesep 'allchnames_' blocknames{i} '.mat']);
        chnames = allchnames.chnames(channels);
    else
        chnames = [];
    end
    
    %% Loop Over Evnt
    
    for j=1:length(evntnames)
        
        ii = unique_index(j); % find whether this has been done before
        
        % Find current evnt name, tell user which sound is processing
        str = evntnames{j};
        disp(['Processing sound ' num2str(j) ': ' str]);
        
        % If the sound file name hasn't been added to out yet, load audio
        if isempty(out(ii).name)
            
            % Load trial #
            out(ii).trial = j;
            
            %% Load soundname into out structure, read audio
            % out.name will be soundname without '.wav'
            % audioread or wavread will load audio and audiof
            
            if length(str)>2 && strcmpi(str(end-2:end),'wav') % if .wav is included in evntname
                
                soundchar = strsplit(str, {'/','.wav','.WAV'}); % split string at .wav or /
                out(ii).name = soundchar{end-1}; % define out.name as name of sound without .wav or /
                
                % Read the audio file
                if exist('audioread','file')
                    [audio,audiof] = audioread([evnt(j).StimPath filesep evnt(j).name]);
                else
                    [audio,audiof] = wavread([evnt(j).StimPath filesep evnt(j).name]);
                end
                
            else % if .wav is not included in evntname
                
                out(ii).name = str; % outname is the same as evnt.name
                
                % Read the audio file
                if exist('audioread','file')
                    [audio,audiof] = audioread([evnt(j).StimPath filesep evnt(j).name,'.wav']);
                else
                    [audio,audiof] = wavread([evnt(j).StimPath filesep evnt(j).name,'.wav']);
                end
                
            end
            
            % Audio sampling rate corrections
            if audiof ~= 16000
                audio = resample(audio,16000,audiof);
                audiof = 16000;
                warning('audiof resampled to 16k for compatibility with wav2aud');
            end
%             if audiof == 44100
%                 audio = resample(audio,16000,audiof);
%                 audiof = 16000;
%                 warning('audiof is 44100; resampled to 16k for to get rid of empty info in spectrogram');
%             end
            
            %% Load audio info into out
            
            % Add bef/aft silence to audio, load to out
            audiodim = min(size(audio)); audiolength = max(size(audio));
            audio = reshape(audio, [audiolength, audiodim]);  % reshape audio to ensure number of channels are second dimension
            audio = [ zeros(round(befaft(1)*audiof),audiodim); audio; zeros(round(befaft(2)*audiof),audiodim) ]; % add the appropriate silence before and after the audio
            out(ii).sound = audio;
            out(ii).soundf = audiof;
            out(ii).befaft = befaft;
            out(ii).duration = evnt(j).stopTime-evnt(j).startTime+sum(befaft);
            
            
            %% Read data and artifact if this is the first time through this block
            
            currblock = evnt(j).block; % check current block in evnt
            
            if strcmp(currblock,prevblock); % if correct block, do nothing
                
            else % if wrong block, load data
                
                if ~isempty(audchannels)
                    % Read audio htk files, build audrecord matrix
                    audrecord = [];
                    for k = 1:length(audchannels)
                        [audrecord(:,k),audrecordf] = readhtk([datapath filesep 'processed' filesep currblock filesep 'analog' filesep 'a' num2str(audchannels(k)) '.htk']);
                    end
                end
                
                % Read data htk files, build data matrix
                tmp = [];
                for k = 1:length(channels)
                    [tmp(:,k),f] = readhtk([datapath filesep 'processed' filesep currblock filesep char(cond) filesep 'Ch' num2str(channels(k)) '.htk']);
                end
                
                % Resample data if indicated
                if no_resample
                    record = tmp;
                    dataf = f;
                else
                    record = resample(tmp,dataf,round(f));
                end
                clear tmp;
 
                % Take reference if present
                reference = [];
                refdir = [datapath filesep 'processed' filesep currblock filesep 'reference'];                
                if exist(refdir,'dir')
                    numfiles = length(dir([refdir filesep 'Ch*']));
                    for k = 1:numfiles;
                        [reference(:,k),f] = readhtk([refdir filesep 'Ch' num2str(channels(k)) '.htk']);
                    end
                end             
                
                % Find artifact: Laura's version has not evaluated this, may still need some adjustment to work with updated code
                if ~isempty(artifact)
                    tmp = [];
                    for k = artifact
                        [tmp(:,k),f_art] = readhtk([datapath filesep currblock '/artifact/' num2str(k)]);
                    end
                    artif_wave = downsample(tmp,round(f_art/dataf));
                end
            end
            
            % Update block
            prevblock = currblock; % update block
            
            
            %% Load data and data info to out
            
            % Load data info to out; set empty for resp and artifact
            out(ii).dataf = dataf; 
            out(ii).type = datatype;
            out(ii).resp = [];
            
            % Load channelnames if present
            if ~isempty(chnames)
                out(ii).chnames = chnames;
            end
            
            % Load audchannels if present
            if ~isempty(audchannels)
                out(ii).audrecord = [];
                out(ii).audrecordf = audrecordf;
            end
            
            % Load artifact if present
            if ~isempty(artifact)
                out(ii).artifact = [];
            end
            
            % Load reference if present
            if ~isempty(reference)
                out(ii).reference = [];
            end
            
                        
            %% Take and load auditory spectrogram
            
            for k = 1:audiodim % if there are extra audio dimensions, get an audio spectrogram for each channel
                
                switch specflag
                    
                    case 'Auditory' % normal auditory spectrogram
                        
                        tmpaud(:,:,k) = wav2aud(audio(:,k), [1000/dataf 1000/dataf -2 log2(audiof/16000)] )'; % take spectrogram, resampled to match with dataf
                        
                        try
                            out(ii).aud(:,:,k) = tmpaud(:,1:round(out(ii).duration*dataf),k); % load specgram into out, cut to be the same length as resp
                        catch
                            warning('Size mismatch')
                            out(ii).aud(:,:,k) = tmpaud(:,:,k); % if there's a size error, disp warning and load full specgram into out
                        end
                        
                        clear tmpaud;
                        
                    case 'AudNorm' % normalized auditory spectrogram (NOT resampled to match with dataf)
                        out(ii).aud(:,:,k) = AUDwrapper(audio(:,k), audiof, 0, 0);
                        
                    case 'mfcc' % mfcc
                        out(ii).mfcc(:,:,k) = melfcc_wrapper(audio(:,k),audiof);

                end
                
            end         
        
        else
            %% Load data for repeat stimuli
            % if out.name is already populated, this is a second instance of the same stimulus, so concatenate the trials
            out(ii).trial = [out(ii).trial j];
        end
         
        %% Cut and load data 
        
        % Cut data
        datastart = round((evnt(j).startTime - befaft(1)) * dataf); % find starting sample
        newdata = record(datastart + 1 : datastart + round(out(ii).duration * dataf), :)'; % cut from starting sample through duration * dataf later
        out(ii).resp = cat(3, out(ii).resp, newdata); % load data (set up to be concatenated with any previous presentations of the same stimulus, but not currently applicable)
        
        % Load artifact to out if present
        if ~isempty(artifact)
            newartifact=artif_wave(datastart+1:datastart+round(out(ii).duration*dataf),:)';
            out(ii).artifact=cat(3,out(ii).artifact,newartifact); % Load artifact (set up to be concatenated with any previous presentations of the same stimulus, but not currently applicable)
        end
        
        if ~isempty(reference)
            newreference = reference(datastart+1:datastart+round(out(ii).duration*dataf),:)';
            out(ii).reference=cat(3,out(ii).reference,newreference);
        end
        
        if ~isempty(audchannels)
            % Cut recorded audio
            audstart = round((evnt(j).startTime - befaft(1)) * audrecordf); % find starting sample
            newaudrecord = audrecord(audstart + 1 : audstart + round(out(ii).duration * audrecordf), :); % cut from starting sample through duration * dataf later
            out(ii).audrecord = cat(3, out(ii).audrecord, newaudrecord); % load data (set up to be concatenated with any previous presentations of the same stimulus, but not currently applicable)
        end
        
        cleancond = strsplit(cond,'_');
        cleanvarfile = [datapath filesep 'processed' filesep blocknames{i} filesep cleancond{1} filesep 'cleaningvars.mat'];
        if exist(cleanvarfile,'file')
            % Cut recorded audio
            load(cleanvarfile);
            out(ii).params = cleanvars;
        end
        
    end
    
    if normout
        silencefile = [datapath filesep 'processed' filesep blocknames{i} filesep 'silencestats.mat'];
        out = zscore_out(out,silencefile);
    end
    
    % Notify user the block is done
    disp(['Completed ' prevblock]);
    
    
    % Ask about saving
    saveevnt = input('Save out? ','s');
    if ~strcmpi(saveevnt,{'y', 'yes', '1'})
    else
        if normout
            savename = [datapath filesep 'processed' filesep prevblock filesep 'out_' subject '_' prevblock '_' cond '_norm.mat'];
        else
            savename = [datapath filesep 'processed' filesep prevblock filesep 'out_' subject '_' prevblock '_' cond '.mat'];
        end
        save(savename, 'out','-v7.3');
        disp(['Saved to ' savename]);
    end
end

% Notify user code has finished
fprintf('Completed NeuralGenOut \n\n');

end
