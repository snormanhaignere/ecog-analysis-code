function plot_example_spectra(signal, sr, good_channels, fig_file)

% Plots example electrode spectra. Used to check preprocessing steps.
% 
% -- Inputs -- 
%  
% signal: [samples x electrode] signal matrix
% 
% sr: sampling rate of the signal
% 
% good_channels: vector with indices of electrodes determined to be of decent
% quality (typically by examining 60 Hz noise)
% 
% fig_file: file to save results to
% 
% 2016-08-15 - Last Modified, Sam NH

n_channels_to_plot = 6;

figure;
% set(gcf, 'Position', [0 0 400 800]);
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
        [px, f] = fftplot2(signal(:,example_good_channels(i)), sr);
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
% export_fig([fig_file '.pdf'], '-pdf', '-transparent');
set(gcf, 'PaperSize', [6 10]);
set(gcf, 'PaperPosition', [0.25 0.25 5.5 9.5]);
print([fig_file '.pdf'],'-dpdf');
print([fig_file '.png'],'-dpng', '-r100');
close all;

