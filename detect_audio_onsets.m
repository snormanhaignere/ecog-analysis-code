function audio_onsets_smps = detect_audio_onsets(audio_signal, sr, ...
    stim_directory, StimOrder, varargin)

% Iteratively finds audio onsets by detecting peaks in the
% cross-correlation function between the audio channel and the expected
% waveforms. After a peak is found surrounding peaks in the
% cross-correlation function are suppressed, and the process repeats until
% no more peaks above a certain tolerance are found (0.4 of the max peak
% aross all timepoints). These peaks are then strung together based on
% the expected isi of the stimuli.
% 
% 2019-11-11: Created, Sam NH

clear I;
I.thresh = 0.5;
I.minisi = [];
I.maxisi = [];
I.breaks = [];
I.win_to_zero = 2;
I.transform = 'none';
I.plot = true;
I.stimstoplot = [];
I.zoomedplot = false;
I.pause = false;
I = parse_optInputs_keyvalue(varargin, I);

n_stim_onsets = length(StimOrder);

if isscalar(I.thresh)
    thresh = ones(n_stim_onsets, 1)*I.thresh;
else
    assert(length(I.thresh)==n_stim_onsets);
    thresh = I.thresh;
end

if I.plot
    if isempty(I.stimstoplot)
        plot_figs = true(n_stim_onsets, 1);
    else
        plot_figs = false(n_stim_onsets, 1);
        plot_figs(I.stimstoplot) = true;
    end
else
    plot_figs = false(n_stim_onsets, 1);
end

% open figure for dynamic plotting
if I.plot
    figh_global = figure;
    set(figh_global, 'Position', [100 100 1000 300]);
    if I.zoomedplot
        figh_zoomed = figure;
        set(figh_zoomed, 'Position', [100 100 1200 800]);
    end
end

cc_peaks = cell(1, n_stim_onsets);

minisi = nan(n_stim_onsets, 1);
for i = 1:n_stim_onsets
    
    % read in audio waveform
    [wav, wavsr] = audioread([stim_directory '/' StimOrder{i} '.wav']);
    
    % minimum isi is the duration of this audio waveform
    minisi(i) = length(wav)/wavsr;
    
    switch I.transform
        case 'none'
            sig = resample(wav, sr, wavsr);
        otherwise
            error('No matching transform for %s', I.transform);
    end
    
    [cc, lag] = xcorr(audio_signal, sig);
    cc = cc/max(cc);
    
    % get rid of negative lags
    xi = lag>=0;
    lag = lag(xi);
    cc = cc(xi);
    clear xi;
    
    cc_orig = cc;
    while any(cc>thresh(i))
        
        % find a peak and add to trigger onsets
        [~,xi] = max(cc);
        best_lag = lag(xi);
        cc_peaks{i} = [cc_peaks{i}, best_lag];
        
        % set surround to zero
        xi = abs(lag - best_lag)/sr < I.win_to_zero;
        cc(xi) = 0;
        
        % plot
        if plot_figs(i)
            if I.zoomedplot
                figure(figh_global);
            end
            clf(figh_global);
            hold on;
            for j = 1:length(cc_peaks{i})
                stem([1 1]*cc_peaks{i}(j)/sr, [0 1]*thresh(i), 'r-', 'LineWidth', 2);
            end
            plot(lag/sr, cc_orig, 'k-');
            plot(lag/sr, ones(size(lag))*thresh(i), 'r--');
            ylim([-1,1]);
            xlim(lag([1,end])/sr);
            xlabel('Time');
            ylabel('Norm CC');
            title(sprintf('stim %d, %s', i, StimOrder{i}));
            drawnow;
        end
    end
    
    n_sound_onsets = length(cc_peaks{i});
    fprintf('Stim %d, %s: %d onsets found\n', i, StimOrder{i}, n_sound_onsets);
    drawnow;
    
    % plot
    if plot_figs(i) && I.zoomedplot
        figure(figh_zoomed);
        clf(figh_zoomed);
        hold on;
        n_rows = n_sound_onsets;
        n_cols = 2;
        simulated_audio = zeros(size(audio_signal));
        for j = 1:n_sound_onsets
            
            subplot(n_rows, n_cols, 1 + (j-1)*2);
            
            % plot audio signal
            win = [-0.1 0.4];
            xi = (round(win(1)*sr):round(win(2)*sr)) + cc_peaks{i}(j);
            plot((xi-cc_peaks{i}(j))/sr*1000, zscore(audio_signal(xi)));
            
            % plot template
            yi = (1:length(sig)) + cc_peaks{i}(j);
            simulated_audio(yi) = simulated_audio(yi) + sig;
            hold on;
            plot((xi-cc_peaks{i}(j))/sr*1000,zscore(simulated_audio(xi)));
            xlim((xi([1,end])-cc_peaks{i}(j))/sr*1000);
            title(sprintf('ons %d', j));
            xlabel('ms')
            legend('Recorded', 'Predicted');
            
            % plot cross-correlation
            subplot(n_rows, n_cols, 2 + (j-1)*2);
            win = [-0.1 0.1];
            xi = (round(win(1)*sr):round(win(2)*sr)) + cc_peaks{i}(j);
            plot((lag(xi)-cc_peaks{i}(j))/sr*1000,cc_orig(xi));
            xlim((lag(xi([1,end]))-cc_peaks{i}(j))/sr*1000);
            ylim([-1,1]);
            title('CC'); xlabel('ms');
        end
    end
    
    if plot_figs(i) && I.pause
        input('Press enter to continue');
    end
end

%% String the peaks together to form a sequence

% determine minimum ISI
if ~isempty(I.minisi)
    if isscalar(I.minisi)
        minisi = ones(n_stim_onsets-1,1) * I.minisi;
    else
        minisi = I.minisi;
    end
end

% determine maximum ISI
if isempty(I.maxisi)
    maxisi = inf(n_stim_onsets-1,1);
else
    if isscalar(I.maxisi)
        maxisi = ones(n_stim_onsets-1,1) * I.maxisi;
    else
        maxisi = I.maxisi;
    end
end

% allow breaks to be arbitrarily long
if ~isempty(I.breaks)
    maxisi(I.breaks) = inf;
end

% string together the peaks based on these min and max ISI values
audio_onsets_smps = min(cc_peaks{1});
for i = 2:n_stim_onsets
    
    % find next peak
    time_range = audio_onsets_smps(end)/sr + [minisi(i-1), maxisi(i-1)];
    xi = cc_peaks{i}/sr > time_range(1) & cc_peaks{i}/sr < time_range(2);
    if sum(xi)==0
        error('No peak found for stim %d, %s', i, StimOrder{i});
    end
    audio_onsets_smps = cat(2, audio_onsets_smps, min(cc_peaks{i}(xi)));
    clear xi;
    
end
assert(length(audio_onsets_smps)==n_stim_onsets);


