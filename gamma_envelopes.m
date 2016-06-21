function [gamma_env, gamma_env_t] = gamma_envelopes(input_signal, input_sr, env_sr, gamma_frequency_range)

% [gamma_env, gamma_env_t] = gamma_envelopes(input_signal, input_sr, env_sr)
% 
% Very simple code for measuring gamma power envelopes from ECoG signal.
% 
% -- Inputs --
% 
% input_signal: [samples x electrode matrix] with the ECoG signals from each electrode
% 
% input_sr: sampling rate of the ecog signals in Hz
% 
% env_sr: sampling rate to downsample the envelopes to
% 
% gamma_frequency_range: frequency range of the gamma power
% 
% -- Outputs -- 
% 
% gamma_env: the downsampled envelopes
% 
% gamma_env_t: corresponding timestamps for the downsampled envelopes
% 
% -- Example --
% bci_run_file = '/Users/svnh2/Dropbox (MIT)/mcdexp-svnh/ecog-analysis-code/example_dataset.dat';
% [input_signal,~,parameters,~] = load_bcidat(bci_run_file);
% input_sr = parameters.SamplingRate.NumericValue;
% env_sr = 100;
% gamma_frequency_range = [70 140];
% [gamma_env, gamma_env_t] = gamma_envelopes(double(input_signal), input_sr, env_sr, gamma_frequency_range);
% plot(gamma_env_t(1:5e3), gamma_env(1:5e3,74));
% 
% 2016-1-26: Created by Sam Norman-Haignere 

% filter parameters
filt_order = 6;

% filter frequencies outside of the gamma power range
fprintf('Buttworth filtering...\n'); drawnow;
[B,A] = butter(filt_order, gamma_frequency_range/(input_sr/2));
gamma_filt = filtfilt(B ,A, input_signal);

% plot the spectrum before and after filtering using fftplot2
addpath('/Users/svnh2/Dropbox (MIT)/mcdexp-svnh/general-audio-code')
figure;
subplot(1,2,1); fftplot2(double(input_signal(:,1)), input_sr); xlim([20 input_sr/2]);
subplot(1,2,2); fftplot2(double(gamma_filt(:,1)), input_sr); xlim([20 input_sr/2]);

% envelope extraction
fprintf('Measuring envelopes...\n'); drawnow;
gamma_env = resample(abs(hilbert(gamma_filt)), env_sr, input_sr);
gamma_env(gamma_env<0) = 0;
gamma_env_t = (0:length(gamma_env)-1)/env_sr;