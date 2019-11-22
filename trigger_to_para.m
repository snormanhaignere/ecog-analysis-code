function trigger_to_para(exp, subjid, r, varargin)

% Use triggers to create the para file with stimulus onsets.
% 
% Wrapper around detect_trigger_onsets.m
% 
% 2019-11-06: Created, Sam NH
% 
% 2019-11-11: Misc edits, Sam NH
% 
% 2019-11-16: Added capacity to save figures

global root_directory;

I.stimprefix = '';
I.stimsuffix = '';
I.trigtemplatefile = [root_directory '/ecog-analysis-code/laura/Pipeline-2019-08-09/NYUtrig.mat'];
I.tol = 0.4;
I.win_to_zero = 1.5;
I.plot = false;
I.dur = NaN;
I.para_suffix = '';
[I,C] = parse_optInputs_keyvalue(varargin, I);

% directory structure and external repositories
project_directory = [root_directory '/' exp];

% directory to save figures to
figure_directory = [project_directory '/figures/trigger-sync/' subjid '/r' num2str(r)];

% load the trigger signal
trigger_MAT_file = [project_directory '/data/ECoG-trigger/' subjid '/r' num2str(r) '.mat'];
load(trigger_MAT_file, 'trigger_signal', 'sr');
assert(size(trigger_signal,2)==1);

% zero exclusion window
if isvar_in_mfile(trigger_MAT_file, 'excludewin')
    load(trigger_MAT_file, 'excludewin');
    if ~isempty(excludewin)
        trigger_signal(excludewin(1):excludewin(2)) = 0;
    end
end

% load the trigger template
X = load(I.trigtemplatefile);
trig_template = resample(X.trigform, sr, X.fs);

% find the trigger onsets
trig_onsets = detect_trigger_onsets(trigger_signal, trig_template, sr, ...
    'tol', I.tol, 'win_to_zero', I.win_to_zero*length(trig_template)/sr, ...
    'plot', I.plot, 'figdir', figure_directory);

% load stimulus orders
stim_order_file = [project_directory '/data/ECoG-stimorders' ...
    '/' subjid '/r' num2str(r) '.mat'];
load(stim_order_file, 'StimOrder','stim_directory');
if length(StimOrder)~=length(trig_onsets)
    error('Wrong number of triggers found.\nExpected %d, found %d\n', ...
        length(StimOrder), length(trig_onsts));
end

% para file to create
if isempty(I.para_suffix)
    para_file = [project_directory '/data/para/' subjid '/r' num2str(r) '.par'];
else
    para_file = [project_directory '/data/para/' subjid '/r' num2str(r) '_' I.para_suffix '.par'];
end

% write para file
fid = fopen(mkpdir(para_file), 'w');
for i = 1:length(StimOrder)
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
        trig_onsets(i)/sr, NaN, dur, NaN, ...
        stimname);
end
fclose(fid);
