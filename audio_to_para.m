function audio_to_para(exp, subjid, r, varargin)

% Use audio channel to create the para file with stimulus onsets.
% 
% Wrapper around detect_audio_onsets.m
% 
% 2019-11-11: Created, Sam NH

global root_directory;

clear I;
I.stimprefix = '';
I.stimsuffix = '';
I.dur = NaN;
I.thresh = 0.5;
I.minisi = [];
I.maxisi = [];
I.breaks = [];
I.win_to_zero = 2;
I.transform = 'none';
I.plot = false;
I.stimstoplot = [];
I.zoomedplot = false;
I.pause = false;
I.para_suffix = '';
[I,C] = parse_optInputs_keyvalue(varargin, I);

% directory structure and external repositories
project_directory = [root_directory '/' exp];

% load the audio signal
audio_MAT_file = [project_directory '/data/ECoG-audio/' subjid '/r' num2str(r) '.mat'];
load(audio_MAT_file, 'audio_signal', 'sr');
assert(size(audio_signal,2));

% zero exclusion window
if isvar_in_mfile(audio_MAT_file, 'excludewin')
    load(audio_MAT_file, 'excludewin');
    if ~isempty(excludewin)
        audio_signal(excludewin(1):excludewin(2)) = 0;
    end
end

% load stimulus orders
stim_order_file = [project_directory '/data/ECoG-stimorders' ...
    '/' subjid '/r' num2str(r) '.mat'];
load(stim_order_file, 'StimOrder', 'stim_directory');
n_stim_onsets = length(StimOrder);

% find the trigger onsets
audio_onsets = detect_audio_onsets(audio_signal, sr, stim_directory, StimOrder, ...
    'thresh', I.thresh, 'minisi', I.minisi, 'maxisi', I.maxisi, ...
    'breaks', I.breaks, 'win_to_zero', I.win_to_zero, ...
    'transform', I.transform, 'plot', I.plot, 'zoomedplot', I.zoomedplot, ...
    'stimstoplot', I.stimstoplot, 'pause', I.pause);

% para file to create
if isempty(I.para_suffix)
    para_file = [project_directory '/data/para/' subjid '/r' num2str(r) '.par'];
else
    para_file = [project_directory '/data/para/' subjid '/r' num2str(r) '_' I.para_suffix '.par'];
end

% write para file
fid = fopen(mkpdir(para_file), 'w');
for i = 1:n_stim_onsets
    stimname = StimOrder{i};
    if C.stimprefix
        stimname = [I.stimprefix '_' stimname]; %#ok<AGROW>
    end
    if C.stimsuffix
        stimname = [stimname '_' I.stimsuffix]; %#ok<AGROW>
    end
    if isnan(I.dur)
        wavinfo = audioinfo([stim_directory '/' StimOrder{i} '.wav']);
        dur = wavinfo.Duration;
    else
        if isscalar(I.dur)
            dur = I.dur;
        else
            assert(length(I.dur)==length(StimOrder));
            dur = I.dur(i);
        end
    end
    fprintf(fid,'%10.6f%5d%10.6f%5d%90s\n', ...
        audio_onsets(i)/sr, NaN, dur, NaN, ...
        stimname);
end
fclose(fid);
