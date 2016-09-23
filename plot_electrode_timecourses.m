function plot_electrode_timecourses(signal, sr, varargin)

% Plots example electrode timecourses. Used to check preprocessing steps.
%
% -- Inputs --
%
% signal: [samples x electrode] signal matrix
%
% sr: sampling rate of the signal
%
% 2016-08-15 - Last Modified, Sam NH
%
% 2016-09-20 - Modified to plot all electrodes and single-mosaic

n_electrodes = size(signal,2);
I.output_file = '';
I.single_mosaic = false;
I.electrode_numbers = 1:n_electrodes;
I = parse_optInputs_keyvalue(varargin, I);

if I.single_mosaic
    
    if n_electrodes > 10
        error('Can''t plot more than 10 electrodes on a single mosaic');
    end
    
    figure;
    set(gcf, 'Position', [0 0 600 150*n_electrodes]);
    for i = 1:n_electrodes
        for j = 1:2
            if j == 1;
                smps_to_plot = 1:size(signal,1);
                titlestring = sprintf(...
                    ['elec ' num2str(I.electrode_numbers(i)) ...
                    ', full-timecourse']);
            else
                midpoint = round(size(signal,1)/2);
                smps_to_plot = midpoint:midpoint+sr-1;
                clear midpoint;
                titlestring = sprintf(...
                    ['elec ' num2str(I.electrode_numbers(i)) ...
                    ', 1-second']);
            end
            
            % plot
            subplot(n_electrodes, 2, j + (i-1)*2);
            plot(smps_to_plot/sr, ...
                signal(smps_to_plot, i), 'k-');
            
            % axis limits
            xlim([smps_to_plot(1), smps_to_plot(end)]/sr);
            clear smps_to_plot;
            
            % axis labels
            set(gca, 'FontSize', 8);
            xlabel('Time (s)'); ylabel('Amplitude');
            title(titlestring);
        end
    end
    
    if ~isempty(I.output_file)
        set(gcf, 'PaperSize', [8, n_electrodes*10/6]);
        set(gcf, 'PaperPosition', [0.25, 0.25, 7.5, n_electrodes*10/6-0.5]);
        print([I.output_file '.pdf'],'-dpdf');
        print([I.output_file '.png'],'-dpng', '-r100');
    end
    
else
    
    figh = figure;
    for i = 1:n_electrodes
        
        clf(figh);
        set(figh, 'Position', [0 0 600 200]);
        
        for j = 1:2
            if j == 1;
                smps_to_plot = 1:size(signal,1);
                titlestring = sprintf(...
                    ['elec ' num2str(I.electrode_numbers(i)) ...
                    ', full-timecourse']);
            else
                midpoint = round(size(signal,1)/2);
                smps_to_plot = midpoint:midpoint+sr-1;
                clear midpoint;
                titlestring = sprintf(...
                    ['elec ' num2str(I.electrode_numbers(i)) ...
                    ', 1-second']);
            end
            
            subplot(1, 2, j);
            plot(smps_to_plot/sr, ...
                signal(smps_to_plot, i), 'k-');
            
            % a limits
            xlim([smps_to_plot(1), smps_to_plot(end)]/sr);
            clear smps_to_plot;
            
            % axis labels
            set(gca, 'FontSize', 8);
            xlabel('Time (s)'); ylabel('Amplitude');
            title(titlestring);
            
        end
        
        if ~isempty(I.output_file)
            set(figh, 'PaperSize', [8 2]);
            set(figh, 'PaperPosition', [0.25 0.25 7.5 1.5]);
            print([I.output_file '-electrode' num2str(I.electrode_numbers(i)) '.pdf'],'-dpdf');
            print([I.output_file '-electrode' num2str(I.electrode_numbers(i)) '.png'],'-dpng', '-r100');
        end
        
    end
    
    
end

