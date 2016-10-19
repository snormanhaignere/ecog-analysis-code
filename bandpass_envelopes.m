function env = bandpass_envelopes(signal, signal_sr, env_sr, ...
    band_in_Hz, filter_order, figure_directory, varargin)

% Measure envelopes from bandpass filters.
% 
% -- Inputs --
% 
% signal: [time x electrode] signal matrix
%
% signal_sr: sampling rate of the signal
% 
% env_sr: sampling rate to use for the envelopes (typically much lower than the
% signal sampling rate, e.g. 100 Hz vs. 1200 Hz)
% 
% band_in_Hz: 2-dimensional vector of frequency cutoffs
% 
% filter_orders: order of the butterworth filters
% 
% figure_directory: directory to save figures to
% 
% -- Outputs --
% 
% envelopes: [time x bands x electrode] matrix with envelopes
% 
% 2016-1-26: Created by Sam NH
% 
% 2016-09-23: Changes to handling of optional inputs and plotting, Sam NH

%% Setup

n_channels = size(signal,2);

I.good_channels = 1:n_channels;
I.electrode_numbers = 1:n_channels;
% I.plot_all_electrodes = false;
I = parse_optInputs_keyvalue(varargin, I);

% directory to save figures with timecourses/spectra for individual electrodes
% single_electrode_directory = [figure_directory '/all-individual-electrodes'];
% if ~exist(single_electrode_directory, 'dir')
%     mkdir(single_electrode_directory);
% end

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass('N,F3dB1,F3dB2', filter_order, ...
    band_in_Hz(1), band_in_Hz(2), signal_sr);
Hd = design(h,'butter');
[B, A] = sos2tf(Hd.sosMatrix,Hd.scaleValues);

% apply filter
fprintf('Applying filter...\n');
subb = filtfilt(B,A,signal);

% measure envelope
fprintf('Calculating envelopes...\n');
env = abs(hilbert(subb));
    
% resample to desired rate
fprintf('Downsampling...\n');
env = resample(env, env_sr, signal_sr);

% truncate
env(env < 0) = 0;

%% Plot subband timecourses and spectra

% example good channels to plot in a mosaic
n_channels_to_plot = 6;
xi = round(linspace(1,length(I.good_channels),n_channels_to_plot+2));
good_channels_to_plot = I.good_channels(xi(2:end-1));
clear xi;

% mosaic of subband timecourses
plot_electrode_timecourses(subb(:,good_channels_to_plot), signal_sr, ...
    'electrode_numbers', I.electrode_numbers(good_channels_to_plot), ...
    'single_mosaic', true, ...
    'output_file',  [figure_directory '/bpfilt_' ...
    num2str(band_in_Hz(1)) '-' num2str(band_in_Hz(2)) 'Hz_subb_timecourse']);

% % subband timecourses for all electrodes
% if I.plot_all_electrodes
%     plot_electrode_timecourses(subb, signal_sr, ...
%         'electrode_numbers', I.electrode_numbers, ...
%         'output_file',  [single_electrode_directory '/bpfilt_' ...
%         num2str(band_in_Hz(1)) '-' num2str(band_in_Hz(2)) 'Hz_subb_timecourse']);
% end

% mosaic of subband spectra
plot_electrode_spectra(subb(:,good_channels_to_plot), signal_sr, ...
    'electrode_numbers', I.electrode_numbers(good_channels_to_plot), ...
    'single_mosaic', true, ...
    'output_file', [figure_directory '/bpfilt_' ...
    num2str(band_in_Hz(1)) '-' num2str(band_in_Hz(2)) 'Hz_subb_spectra']);
    
% subband spectra for all electrodes
% if I.plot_all_electrodes
%     plot_electrode_timecourses(subb, signal_sr, ...
%         'electrode_numbers', I.electrode_numbers, ...
%         'output_file',  [single_electrode_directory '/bpfilt_' ...
%         num2str(band_in_Hz(1)) '-' num2str(band_in_Hz(2)) 'Hz_subb_spectra']);
% end

%% Plot envelope timecourses and spectra

% mosaic of envelope timecourses
plot_electrode_timecourses(env(:,good_channels_to_plot), env_sr, ...
    'electrode_numbers', I.electrode_numbers(good_channels_to_plot), ...
    'single_mosaic', true, ...
    'output_file',  [figure_directory '/bpfilt_' ...
    num2str(band_in_Hz(1)) '-' num2str(band_in_Hz(2)) 'Hz_env_timecourse']);

% envelope timecourses for all electrodes
% if I.plot_all_electrodes
%     plot_electrode_timecourses(env, signal_sr, ...
%         'electrode_numbers', I.electrode_numbers, ...
%         'output_file',  [single_electrode_directory '/bpfilt_' ...
%         num2str(band_in_Hz(1)) '-' num2str(band_in_Hz(2)) 'Hz_env_timecourse']);
% end

