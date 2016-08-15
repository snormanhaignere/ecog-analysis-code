function plot_example_timecourses(signal, sr, good_channels, fig_file)

% Plots example electrode timecourses. Used to check preprocessing steps.
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
set(gcf, 'Position', [0 0 600 800]);
example_good_channels = round(linspace(1,length(good_channels),n_channels_to_plot+2));
example_good_channels = example_good_channels(2:end-1);
for i = 1:n_channels_to_plot
    for j = 1:2
        if j == 1;
            smps_to_plot = 1:size(signal,1);
            titlestring = sprintf(...
                ['elec ' num2str(example_good_channels(i)) ...
                ', full-timecourse']);
        else
            midpoint = round(size(signal,1)/2);
            smps_to_plot = midpoint:midpoint+sr-1;
            clear midpoint;
            titlestring = sprintf(...
                ['elec ' num2str(example_good_channels(i)) ...
                ', 1-second']);
        end
        
        % plot
        subplot(n_channels_to_plot, 2, j + (i-1)*2);
        plot(smps_to_plot/sr, ...
            signal(smps_to_plot, example_good_channels(i)), 'k-');
        
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
% export_fig([fig_file '.pdf'], '-pdf', '-transparent');
set(gcf, 'PaperSize', [8 10]);
set(gcf, 'PaperPosition', [0.25 0.25 7.5 9.5]);
print([fig_file '.pdf'],'-dpdf');
print([fig_file '.png'],'-dpng', '-r100');