function [out] = NeuralGenOut_Laura(datapath, cond, subject, blocks, channels, dataf, befaft, specflag, datatype, artifact)
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
    channels = getchannelinds([datapath filesep 'processed' filesep blocknames{1} filesep 'htkraw'],'htk'); % find channel inds
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
if ~exist('artifact','var')
    artifact=[];
end

%% Loop over blocks

for i = 1:length(blocknames)
    
    disp(['Processing ' blocknames{i} '...']);
    
    % Reset evnt, out, data; load the block's evnt
    evnt = []; out = []; record = [];
    load([datapath filesep 'Processed' filesep blocknames{i} filesep 'evnt_' subject, '_' blocknames{i}, '.mat']);
    
    % Load evntnames, make out.name appropriate length
    evntnames = {evnt.name}';
    out = struct('name',cell(1,length(evntnames)));
    
    prevblock='  '; % previous block; keeps track of whether the block has changed and new data should be loaded; starts blank
    loadload; close; % loads audio colormap, filterbank, params, gets rid of resulting figure
    
    %% Loop Over Evnt
    
    for j=1:length(evntnames)
        
        % If the sound file name hasn't been added to out yet, load audio
        if isempty(out(j).name)
            
            % Load trial #
            out(j).trial = j;
            
            % Find current evnt name, tell user which sound is processing
            str = evntnames{j};
            disp(['Processing sound ' num2str(j) ': ' str]);
            
            %% Load soundname into out structure, read audio
            % out.name will be soundname without '.wav'
            % audioread or wavread will load audio and audiof
            
            if strcmpi(str(end-2:end),'wav') % if .wav is included in evntname
                
                soundchar = strsplit(str, {'/','.wav','.WAV'}); % split string at .wav or /
                out(j).name = soundchar{end-1}; % define out.name as name of sound without .wav or /
                
                % Read the audio file
                if exist('audioread')
                    [audio,audiof] = audioread([evnt(j).StimPath filesep evnt(j).name]);
                else
                    [audio,audiof] = wavread([evnt(j).StimPath filesep evnt(j).name]);
                end
                
            else % if .wav is not included in evntname
                
                out(j).name = str; % outname is the same as evnt.name
                
                % Read the audio file
                if exist('audioread')
                    [audio,audiof] = audioread([evnt(j).StimPath filesep evnt(j).name,'.wav']);
                else
                    [audio,audiof] = wavread([evnt(j).StimPath filesep evnt(j).name,'.wav']);
                end
                
            end
            
            %% Load audio info into out
            
            % Add bef/aft silence to audio, load to out
            audio = audio(:,1); % remove any extra dimensions from audio
            audio = [ zeros(befaft(1)*audiof,1); audio; zeros(befaft(2)*audiof,1) ]; % add the appropriate silence before and after the audio
            out(j).sound = audio;
            out(j).soundf = audiof;
            out(j).befaft = befaft;
            out(j).duration = evnt(j).stopTime-evnt(j).startTime+sum(befaft);
            
            
            %% Read data and artifact if this is the first time through this block
            
            currblock = evnt(j).block; % check current block in evnt
            
            if strcmp(currblock,prevblock); % if correct block, do nothing
                
            else % if wrong block, load data
                
                % Read data htk files, build data matrix
                tmp = [];
                for k = 1:length(channels)
                    [tmp,f] = readhtk([datapath filesep 'processed' filesep currblock filesep char(cond) filesep 'Ch' num2str(channels(k)) '.htk']);
                    for cnt2 = 1:size(tmp,1)
                        if cnt2==1,
                            fieldname = 'resp';
                        else
                            fieldname = ['resp' num2str(cnt2)];
                        end
                        tmp2 = tmp(cnt2,:);
                    end
                end
                % Resample data if indicated
                if no_resample
                    record = tmp;
                    dataf = f;
                else
                    record = resample(tmp,dataf,round(f));
                end
                
                % Find artifact: Laura's version has not evaluated this, may still need some adjustment to work with updated code
                tmp = [];
                for k = artifact
                    [tmp(:,k),f_art] = readhtk([datapath filesep currblock '/artifact/' num2str(k)]);
                end
                if ~isempty(artifact)
                    warning('Laura''s editing note: artifact code not yet checked');
                    artif_wave = downsample(tmp,round(f_art/dataf));
                end
            end
            
            % Update block
            prevblock = currblock; % update block
            
            
            %% Load data and data info to out
            
            % Load data info to out; set empty for resp and artifact
            out(j).dataf = dataf; 
            out(j).type = datatype;
            out(j).resp = [];
            out(j).artifact = []; % Could be moved to 
            
            % Load channel names if they're in evnt
            if isfield(evnt,'channelnames') % loads the channelnames if current evnt isn't written to contain them
                out(j).channelnames=evnt(j).channelnames;
            end
            
            % Cut data
            datastart = round((evnt(j).startTime - befaft(1)) * dataf); % find starting sample
            newdata = record(datastart + 1 : datastart + round(out(j).duration * dataf), :)'; % cut from starting sample through duration * dataf later
            out(j).resp = cat(3, out(j).resp, newdata); % load data (set up to be concatenated with any previous presentations of the same stimulus, but not currently applicable)
            
            % Load artifact to out if present
            if ~isempty(artifact)
                newartifact=artif_wave(datastart+1:datastart+round(out(j).duration*dataf),:)';
                out(j).artifact=cat(3,out(j).artifact,newartifact); % Load artifact (set up to be concatenated with any previous presentations of the same stimulus, but not currently applicable)
            end
            
                        
            %% Take and load auditory spectrogram
            
            switch specflag
                
                case 'Auditory' % normal auditory spectrogram
                    tmpaud = wav2aud(audio, [1000/dataf 1000/dataf -2 log2(audiof/16000)] )'; % take spectrogram, resampled to match with dataf
                    
                    try
                        out(j).aud=tmpaud(:,1:round(out(j).duration*dataf)); % load specgram into out, cut to be the same length as resp
                    catch
                        warning('Size mismatch')
                        out(j).aud=tmpaud; % if there's a size error, disp warning and load full specgram into out
                    end
                    
                case 'AudNorm' % normalized auditory spectrogram (NOT resampled to match with dataf)
                    out(j).aud = AUDwrapper(audio, audiof, 0, 0);
                        
                case 'mfcc' % mfcc
                    out(j).mfcc=melfcc_wrapper(audio,audiof);
            end                       
            
        
        else
            %% Load data for repeat stimuli
            % if out.name is already populated, this is a second instance
            % of the same stimulus, so we need to concatenate the trials;
            % however I don't think I'm set up for this so need to go back
            % and make it an option
            warning('Laura''s editing note: double-check code for cases where a stimulus is presented multiple times; may not be supported currently');
            out(j).trial = [out(j).trial j];
        end
        
    end
    
    % Notify user the block is done
    disp(['Completed ' prevblock]);
    
    % Ask about saving
    saveevnt = input('Save out? ','s');
    if ~strcmpi(saveevnt,{'y', 'yes', '1'})
    else
        save([datapath filesep 'Processed' filesep prevblock filesep 'out_' subject '_' prevblock '_' cond '.mat'], 'out','-v7.3')
        disp(['Saved to ' datapath filesep 'Processed' filesep prevblock filesep 'out_' subject '_' prevblock '_' cond '.mat']);
    end
end

% Notify user code has finished
fprintf('Completed NeuralGenOut \n\n');

end
