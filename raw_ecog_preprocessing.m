function [signal, good_channels, noise60Hz_rms] = ...
    raw_ecog_preprocessing(signal, sr, figure_directory, varargin)

% Preprocessing steps applied to the raw ecog data. Steps include
%
% (1) Determines which channels are good/bad based on deviations in 60 Hz line noise
% across electrodes.
%
% (2) High-pass filters the signal with a very low cutoff (0.5
% Hz).
%
% (3) Remove the common average reference of good channels.
%
% (4) Notch filter the signal to remove 60 Hz noise.
%
% Plots timecourses and spectra after each step.
%
% -- Input parameters --
%
% signal: [time x electrode] signal matrix
%
% sr: signal sampling rate in Hz
%
% figure_directory: directory to save plots to
%
% figure_fname_prefix: file name prefix for all plots
%
% 2016-07-19, Created by Sam NH
%
% 2016-08-12, Generalized to work with an arbitrary input signal and coded a
% separate wrapper to handle data input/output
% 
% 2019-01-21: Changed parameter handling

% number of total channels
n_channels = size(signal,2);

% optional variables / inputs
I.steps = {'60Hz', 'car', 'notch'};
I.bw60Hz = 0.6;
I.frac = 0.5; % see channel_selection_from_60Hz_noise.m
I.min_nchannels = 10;
I.thresh60Hz = 5;
I.hpcutoff = 0.5;
I.hporder = 4;
I.notchbw = 1;
I.notchfreqs = [60, 120, 180, 240];
I.electrode_numbers = 1:n_channels;
I.plot_all_electrodes = false;
I.exclude_from_car = [];
I.array_inds = [];
I.chnames = {};
I = parse_optInputs_keyvalue(varargin, I);

% directory to save figures with timecourses for individual electrodes
single_electrode_directory = [figure_directory '/all-individual-electrodes'];
if ~exist(single_electrode_directory, 'dir')
    mkdir(single_electrode_directory);
end

% detect 60 Hz noise
if any(strcmp('60Hz', I.steps))
    fprintf('Detecting good/bad channels with 60 Hz noise...\n');
    [good_channels, noise60Hz_rms] = ...
        channel_selection_from_60Hz_noise(signal, sr, ...
        [figure_directory '/60Hz_noise'], 'bw', I.bw60Hz, ...
        'frac', I.frac, 'thresh', I.thresh60Hz, 'min_nchannels', I.min_nchannels,...
        'chnames', I.chnames);
    I.steps = setdiff(I.steps, '60Hz');
elseif any(strcmp('array_60Hz', I.steps))
    assert(~any(strcmp('60Hz', I.steps)));
    fprintf('Detecting good/bad channels using 60 Hz noise within arrays...\n');
    [good_channels, noise60Hz_rms] = ...
        channel_selection_from_60Hz_noise(signal, sr, ...
        [figure_directory '/60Hz_noise_array'], 'bw', I.bw60Hz, ...
        'frac', I.frac, 'thresh', I.thresh60Hz, 'array_inds', I.array_inds, ...
        'min_nchannels', I.min_nchannels, 'chnames', I.chnames);
    I.steps = setdiff(I.steps, 'array_60Hz');
else
    good_channels = 1:n_channels;
    noise60Hz_rms = [];
end

% example good channels to plot in a mosaic
n_channels_to_plot = 6;
xi = round(linspace(1,length(good_channels),n_channels_to_plot+2));
good_channels_to_plot = good_channels(xi(2:end-1));
clear xi;

% plot electrode timecourses in a mosaic, and exhaustively for all electrodes
plot_electrode_timecourses(signal(:,good_channels_to_plot), sr, ...
    'electrode_numbers', I.electrode_numbers(good_channels_to_plot), ...
    'single_mosaic', true, ...
    'output_file', [figure_directory '/timecourse_0-raw']);
if I.plot_all_electrodes
    plot_electrode_timecourses(signal, sr, ...
        'output_file', [single_electrode_directory '/timecourse_0-raw']);
end

% plot electrode spectra in a mosaic, and exhaustively for all electrodes
plot_electrode_spectra(signal(:,good_channels_to_plot), sr, ...
    'electrode_numbers', I.electrode_numbers(good_channels_to_plot), ...
    'single_mosaic', true, ...
    'output_file', [figure_directory '/spectrum_0-raw']);
if I.plot_all_electrodes
    plot_electrode_spectra(signal, sr, ...
        'output_file', [single_electrode_directory '/spectrum_0-raw']);
end
n_steps = length(I.steps);
for i = 1:n_steps
    
    switch I.steps{i}
        
        case 'highpass'
            
            fprintf('High-pass filtering the signal...\n');
            signal = hp_filt(signal, sr, 'order', I.hporder, 'cutoffs', I.hpcutoff);
            
        case 'car'
            
            fprintf('Subtracting the average timecourse of good channels...\n');
            common_average = mean(signal(:,setdiff(good_channels, I.exclude_from_car)),2);
            signal = signal - common_average * ones(1,n_channels);
            
        case 'array_car'
            
            fprintf('Subtracting the average timecourse of good channels within each array...\n');
            assert(length(I.array_inds)==size(signal,2));
            unique_inds = unique(I.array_inds);
            for j = 1:length(unique_inds)
                array_chan = find(unique_inds(j)==I.array_inds);
                chan_to_average = setdiff(intersect(array_chan, good_channels), I.exclude_from_car);
                assert(~isempty(chan_to_average));
                common_average = mean(signal(:,chan_to_average),2);
                signal(:,array_chan) = bsxfun(@minus, signal(:,array_chan), common_average);
            end
            clear unique_inds array_chan common_average chan_to_average;
            
        case 'notch'
            
            fprintf('Notch filtering line noise...\n');
            signal = notch_filt(signal, sr, 'bw', I.notchbw, 'freqs', I.notchfreqs);
            
        otherwise
            
            error('No matching step for %s\n', I.steps{i});
            
    end
    
    % plot electrode timecourses in a mosaic, and exhaustively for all electrodes
    plot_electrode_timecourses(signal(:,good_channels_to_plot), sr, ...
        'electrode_numbers', I.electrode_numbers(good_channels_to_plot), ...
        'single_mosaic', true, ...
        'output_file', [figure_directory '/timecourse_' num2str(i) '-' I.steps{i}]);
    if I.plot_all_electrodes
        plot_electrode_timecourses(signal, sr, ...
            'output_file', [single_electrode_directory '/timecourse_' num2str(i) '-' I.steps{i}]);
    end
    
    % plot electrode spectra in a mosaic, and exhaustively for all electrodes
    plot_electrode_spectra(signal(:,good_channels_to_plot), sr, ...
        'electrode_numbers', I.electrode_numbers(good_channels_to_plot), ...
        'single_mosaic', true, ...
        'output_file', [figure_directory '/spectrum_' num2str(i) '-' I.steps{i}]);
    if I.plot_all_electrodes
        plot_electrode_spectra(signal, sr, ...
            'output_file', [single_electrode_directory '/spectrum_' num2str(i) '-' I.steps{i}]);
    end
end





