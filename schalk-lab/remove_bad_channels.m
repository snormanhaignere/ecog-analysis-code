% Author: Joshua Marple (marpljos@gmail.com), July 2014
%
% This code is part of a project that is meant to provide initial data for
% Dr. Schalk's function through bias modulation hypothesis.
%
% Description:
%   This function tells you the indices of the good channels of your data
%   set. It finds them using an iir peak filter. It also returns to you the
%   noise level of each channel, 
%
% Input:
%   signal - the signal you wish to analyze
%   sampling_rate - the sampling rate of the signal
%
% Output:
%   good_chs - a list of the good channels
%   noise_level - a vector containing the line noise of the good channels
%
% Example:
%   [good_chs, noise] = remove_bad_channels(signal, 1200);

function [good_chs, noise_level] = remove_bad_channels(signal, sampling_rate)

    % iir peak filter
    num_chs = size(signal, 2);
    [b,a] = iirpeak(60/(sampling_rate/2), 0.001);
    noise_level = mean(sqrt(filter(b,a, signal).^2),1);
    noise_mean = mean(noise_level(:, 1:num_chs));
    std_dev = std(noise_level(:, 1:num_chs));
    good_chs = [];
    
%     for idx = 1:num_chs
%         if (noise_level(idx) > noise_mean - std_dev) && (noise_level(idx) < noise_mean + std_dev)
%             good_chs = [good_chs idx]; %#ok<AGROW>
%         end
%     end
    
    for idx = 1:num_chs
        if (noise_level(idx) > median(noise_level) - 10*mad(noise_level, 1)) && (noise_level(idx) < median(noise_level) + 10*mad(noise_level, 1))
            good_chs = [good_chs idx]; %#ok<AGROW>
        end
    end
    
    noise_level = noise_level;
end