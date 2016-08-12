function raw_ecog_preprocessing(schalk_subjid, varargin)

% Preprocessing scripts applied to the raw ecog data.
%
% (1)Determines which channels are good/bad based on deviations in 60 Hz line noise
% across electrodes.
%
% (2) High-pass filters the signal with a very low cutoff (0.5
% Hz).
%
% (3) Removes the common average reference of the good channels.
%
% (4) Notch filters the signal to remove 60 Hz noise.
%
% Plots timecourses and spectra after each step.

% 2016-07-19, Created by Sam NH

% general-purpose ecog analysis code
addpath([root_directory '/ecog-analysis-code']);
addpath([root_directory '/general-experiment-code']);
addpath([root_directory '/schalk-ecog-code']);

% directory for this project
project_directory = [root_directory '/naturalsound-ecog'];

% directory with the data for this experiment
data_directory = [project_directory '/data/ECoG/' schalk_subjid '/'];

% directory to save results to
analysis_directory = [project_directory '/analysis' ...
    '/raw-ecog-preprocessing/' schalk_subjid];
if ~exist(analysis_directory, 'dir');
    mkdir(analysis_directory);
end

% directory to save figures to
figure_directory = strrep(analysis_directory,'analysis','figures');
if ~exist(figure_directory, 'dir');
    mkdir(figure_directory);
end

% parse run numbers from files in data directory
files_in_data_directory = mydir(data_directory);
runs = [];
for i = 1:length(files_in_data_directory)
    runstr = regexp(files_in_data_directory{i}, 'r(\d)+', 'match');
    if ~isempty(runstr)
        runs = [runs, str2double(runstr{1}(2:end))]; %#ok<AGROW>
    end
end
n_runs = length(runs);

%% Load data

for q = 5%:n_runs;
        
    matfile = [analysis_directory '/r' num2str(runs(q)) '.mat'];
    if exist(matfile, 'file') && ~optInputs(varargin, 'overwrite')
        continue
    end
    
    bci_run_file = [data_directory '/r' num2str(runs(q)) '.dat'];
    
    % load the raw data and parameters
    % convert to double precision
    fprintf('Loading signal...\n'); drawnow;
    [raw_signal, ~, params, ~] = load_bcidat(bci_run_file);
    raw_signal = double(raw_signal);
    ecog_sr = params.SamplingRate.NumericValue;
    n_channels = size(raw_signal,2);
    
    %% Select good channels based on 60 Hz noise
    
    % 60-Hz noise stats
    fprintf('Measuring 60 Hz noise...\n'); drawnow;
    [b,a] = iirpeak(60/(ecog_sr/2), 0.001);
    noise60Hz_rms = mean(sqrt(filter(b, a, raw_signal).^2),1);
    noise60Hz_deviation = noise60Hz_rms - median(noise60Hz_rms);
    noise60Hz_thresh = 10*mad(noise60Hz_deviation, 1);
    
    % plot noise
    figure;
    bar(1:n_channels,noise60Hz_deviation);
    hold on;
    plot([1,n_channels], noise60Hz_thresh*[1 1], 'r--');
    plot([1,n_channels], -noise60Hz_thresh*[1 1], 'r--');
    xlabel('Electrodes'); ylabel('RMS Power (Dev. from Median)')
    box off;
    export_fig([figure_directory '/r' num2str(runs(q)) '_60Hz_noise.pdf'],...
        '-pdf', '-transparent');
    close all;
        
    % select good channels
    good_channels = find(abs(noise60Hz_deviation) < noise60Hz_thresh);
    
    % plot timecourse and spectrum of a few good channels
    plot_example_timecourses(raw_signal, ecog_sr, good_channels, ...
        [figure_directory '/r' num2str(runs(q)) '_timecourse_0-raw_sig']);
    plot_example_spectra(raw_signal, ecog_sr, good_channels, ...
        [figure_directory '/r' num2str(runs(q)) '_spectrum_0-raw_sig']);
    
    %% HP filter
    
    % high-pass filter the signal
    hp_filter_order = 4;
    hp_cutoff_in_Hz = 0.5;
    Hd = design(fdesign.highpass(...
        'N,F3dB', hp_filter_order, hp_cutoff_in_Hz, ecog_sr), ...
        'butter');
    [b, a] = sos2tf(Hd.sosMatrix,Hd.scaleValues);
    hp_signal = filtfilt(b,a,raw_signal);
    clear b a Hd;
    
    % plot timecourse and spectrum of high-passed signal
    plot_example_timecourses(hp_signal, ecog_sr, good_channels, ...
        [figure_directory '/r' num2str(runs(q)) '_timecourse_1-after-hpfilter']);
    plot_example_spectra(hp_signal, ecog_sr, good_channels, ...
        [figure_directory '/r' num2str(runs(q)) '_spectrum_1-after-hpfilter']);
    
    %% Common-average reference
    
    % remove common average reference
    fprintf('Removing common average reference...\n');
    common_average = mean(hp_signal(:,good_channels),2);
    car_signal = hp_signal - common_average * ones(1,n_channels);
    
    % plot timecourse and spectrum of common average referenced signal
    plot_example_timecourses(car_signal, ecog_sr, good_channels, ...
        [figure_directory '/r' num2str(runs(q)) '_timecourse_2-after-comm-av-ref']);
    plot_example_spectra(car_signal, ecog_sr, good_channels, ...
        [figure_directory '/r' num2str(runs(q)) '_spectrum_2-after-comm-av-ref']);
    
    %% Notch filtering
    
    fprintf('Notch filtering 60 Hz line noise and harmonics...\n');
    
    % notch filter parameters
    n_harmonics = 4;
    bw = 1/(ecog_sr/2);
    notch.fcenter = (1:n_harmonics)*60;
    notch.bw = ones(length(notch.fcenter),1).*bw;
    notch.wo = notch.fcenter./(ecog_sr/2);
    notch.b = cell(1,n_harmonics); notch.a = cell(1,n_harmonics);
    for i = 1:n_harmonics
        [notch.b{i},notch.a{i}] = iirnotch(notch.wo(i),notch.bw(i));
    end
    
    % fvtool(notch.b{2}, notch.a{2}, 'Fs', ecog_sr);
    % figure;
    % [h,t] = impz(notch.b{1},notch.a{1},[],ecog_sr);
    % plot(t,h);
    % xlim([0 1])
    
    % apply notch filter
    notch_filt_signal = car_signal;
    for i = 1:n_channels,
        for j = 1:n_harmonics,
            notch_filt_signal(:,i) = filtfilt(...
                notch.b{j},notch.a{j},notch_filt_signal(:,i));
        end
    end
    
    % plot timecourse and spectrum of common average referenced signal
    plot_example_timecourses(notch_filt_signal, ecog_sr, good_channels, ...
        [figure_directory '/r' num2str(runs(q)) '_timecourse_3-after-notch-filter']);
    plot_example_spectra(notch_filt_signal, ecog_sr, good_channels, ...
        [figure_directory '/r' num2str(runs(q)) '_spectrum_3-after-notch-filter']);
    
    % save
    preproc_signal = notch_filt_signal;
    save(matfile, 'preproc_signal', 'ecog_sr', 'good_channels');
    
end

function plot_example_timecourses(sig, sr, good_channels, fig_file)

n_channels_to_plot = 6;

figure;
set(gcf, 'Position', [0 0 600 800]);
example_good_channels = round(linspace(1,length(good_channels),n_channels_to_plot+2));
example_good_channels = example_good_channels(2:end-1);
for i = 1:n_channels_to_plot
    for j = 1:2
        if j == 1;
            smps_to_plot = 1:size(sig,1);
            titlestring = sprintf(...
                ['elec ' num2str(example_good_channels(i)) ...
                ', full-timecourse']);
        else
            midpoint = round(size(sig,1)/2);
            smps_to_plot = midpoint:midpoint+sr-1;
            clear midpoint;
            titlestring = sprintf(...
                ['elec ' num2str(example_good_channels(i)) ...
                ', 1-second']);
        end
        
        % plot
        subplot(n_channels_to_plot, 2, j + (i-1)*2);
        plot(smps_to_plot/sr, ...
            sig(smps_to_plot, example_good_channels(i)), 'k-');
        
        % a limits
        xlim([smps_to_plot(1), smps_to_plot(end)]/sr);
        clear smps_to_plot;
        
        % axis labels
        set(gca, 'FontSize', 8);
        xlabel('Time (s)'); ylabel('Amplitude');
        title(titlestring);
    end
end
box off;
export_fig([fig_file '.pdf'], '-pdf', '-transparent');
close all;

function plot_example_spectra(sig, sr, good_channels, fig_file)

n_channels_to_plot = 6;

figure;
set(gcf, 'Position', [0 0 400 800]);
example_good_channels = round(linspace(1,length(good_channels),n_channels_to_plot+2));
example_good_channels = example_good_channels(2:end-1);
for i = 1:n_channels_to_plot
    for j = 1:2
        if j == 1;
            freq_to_plot = [0 sr/2];
            xticks = [0.01 0.1 1 10 100 1000];
            titlestring = sprintf(...
                ['elec ' num2str(example_good_channels(i)) ...
                ', full-spec']);
        else
            freq_to_plot = [50 200];
            xticks = [60 120 180];
            titlestring = sprintf(...
                ['elec ' num2str(example_good_channels(i)) ...
                ', higher-freq']);
        end
        
        % plot
        subplot(n_channels_to_plot, 2, j + (i-1)*2);
        [px, f] = fftplot2(sig(:,example_good_channels(i)), sr);
        semilogx(f, 10*log10(px));
        xlim(freq_to_plot);
        
        % axis labels
        set(gca, 'FontSize', 8);
        set(gca, 'XTick', xticks);
        xlabel('Freq (Hz)'); ylabel('Power (dB)');
        title(titlestring);
    end
end

box off;
export_fig([fig_file '.pdf'], '-pdf', '-transparent');
close all;







