function [signal, good_channels, noise60Hz_rms] = raw_ecog_preprocessing(...
    signal, sr, P, figure_directory, figure_fname_prefix, varargin)

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
% P: parameter structure, see preprocessing_parameters.m
% 
% figure_directory: directory to save plots to
% 
% figure_fname_prefix: file name prefix for all plots
%
% 2016-07-19, Created by Sam NH
%
% 2016-08-12, Generalized to work with an arbitrary input signal and coded a
% separate wrapper to handle data input/output

% by default assume all channels are good
n_channels = size(signal,2);

% preprocessing steps
steps = {'60Hz', 'highpass', 'car', 'notch'};
if optInputs(varargin, 'steps');
    steps = varargin{optInputs(varargin, 'steps')+1};
end
n_steps = length(steps);

% detect 60 Hz noise
if any(strcmp('60Hz', steps))
    fprintf('Detecting good/bad channels with 60 Hz noise...\n');
    [good_channels, noise60Hz_rms] = ...
        channel_selection_from_60Hz_noise(signal, P, sr, ...
        [figure_directory '/' figure_fname_prefix '_60Hz_noise']);
    steps = setdiff(steps, '60Hz');
    n_steps = length(steps);
else
    good_channels = 1:n_channels;
    noise60Hz_rms = [];
end

% plot example timecourses and spectra
% from a few good channels before any preprocessing
plot_electrode_timecourses(signal, sr, good_channels, ...
    [figure_directory '/' figure_fname_prefix ...
    '_timecourse_0-raw']);

plot_electrode_spectra(signal, sr, good_channels, ...
    [figure_directory '/' figure_fname_prefix ...
    '_spectrum_0-raw']);

for i = 1:n_steps
    
    switch steps{i}
            
        case 'highpass'
            
            fprintf('High-pass filtering the signal...\n');
            signal = hp_filt(signal, P, sr);
            
        case 'car'
            
            fprintf('Subtracting the average timecourse of good channels...\n');
            common_average = mean(signal(:,good_channels),2);
            signal = signal - common_average * ones(1,n_channels);

        case 'notch'
            
            fprintf('Notch filtering line noise...\n');
            signal = notch_filt(signal, P, sr);
            
        otherwise
            
            error('No matching step for %s\n', steps{i});
            
    end
    
    % plot timecourses and spectra after each preprocessing step
    plot_example_timecourses(signal, sr, good_channels, ...
        [figure_directory '/' figure_fname_prefix ...
        '_timecourse_' num2str(i) '-' steps{i}]);
    
    plot_example_spectra(signal, sr, good_channels, ...
        [figure_directory '/' figure_fname_prefix ...
        '_spectrum_' num2str(i) '-' steps{i}]);
    
end

close all;






