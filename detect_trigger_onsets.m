function trig_onsets = detect_trigger_onsets(trigger_signal, trig_template, sr, varargin)

% Iteratively finds trigger onsets by detecting peaks in the
% cross-correlation function between the trigger channel and a known
% template. After a peak is found surrounding peaks in the
% cross-correlation function are suppressed, and the process repeats until
% no more peaks above a certain tolerance are found (0.4 of the max peak
% aross all timepoints).
% 
% 2019-06-11: Created, Sam NH

clear I;
I.tol = 0.4;
I.win_to_zero = length(trig_template)/sr*1.5;
I.plot = true;
I = parse_optInputs_keyvalue(varargin, I);

% delays = [1, 3000, 6000]-1;
% fake_trigger_signal = zeros(size(trigger_signal));
% for i = 1:length(delays)
%     xi = (1:length(trig_template))+delays(i);
%     fake_trigger_signal(xi) = trig_template;
% end

% normalized cross-correlation
[cc, lag] = xcorr(trigger_signal, trig_template);
cc = cc/max(cc);

% get rid of negative lags
xi = lag>=0;
lag = lag(xi);
cc = cc(xi);
clear xi;

% open figure for dynamic plotting
if I.plot
    figh = figure;
    set(figh, 'Position', [100 100 1000 300]);
end

% iteratively find triggers
trig_onsets = [];
while max(cc) > I.tol
    
    % find a peak and add to trigger onsets
    [~,xi] = max(cc);
    best_lag = lag(xi);
    trig_onsets = [trig_onsets, best_lag]; %#ok<AGROW>
    
    % set surround to zero
    xi = abs(lag - best_lag)/sr < I.win_to_zero;
    cc(xi) = 0;
    
    % plot
    if I.plot
        clf(figh);
        plot(lag/sr, cc);
        hold on;
        for i = 1:length(trig_onsets)
            plot([1 1]*trig_onsets(i)/sr, [0 1], 'r-', 'LineWidth', 2);
        end
        ylim([0,1]);
        xlabel('Time');
        ylabel('Norm CC');
        drawnow;
    end
    
end

% print number of triggers found
trig_onsets = sort(trig_onsets);
fprintf('Found %d triggers\n', length(trig_onsets));
drawnow;

if I.plot
    
    % plot simulated trigger signal
    n_triggers = length(trig_onsets);
    simulated_trigger_signal = zeros(size(trigger_signal));
    for i = 1:n_triggers
        xi = (1:length(trig_template)) + trig_onsets(i);
        simulated_trigger_signal(xi) = trig_template;
    end
    simulated_trigger_signal = max(trigger_signal)*simulated_trigger_signal/max(simulated_trigger_signal);
    
    figh = figure;
    set(figh, 'Position', [100 100 1000 300]);
    plot((0:length(trigger_signal)-1)/sr, [trigger_signal, simulated_trigger_signal]);
    xlabel('Time');
    legend('Actual', 'Predicted', 'Location', 'EastOutside');
    
    figh = figure;
    set(figh, 'Position', [100 100 1000 600]);
    trigs_to_plot = round(linspace(1, n_triggers, 4+2));
    trigs_to_plot = trigs_to_plot(2:5);
    for i = 1:length(trigs_to_plot)
        subplot(2,2,i);
        xi = (1-sr:length(trig_template)+sr) + trig_onsets(trigs_to_plot(i));
        plot((xi-1)/sr, [trigger_signal(xi), simulated_trigger_signal(xi)]);
        xlabel('Time');
        title(sprintf('trigger %d', trigs_to_plot(i)));
        legend('Actual', 'Predicted', 'Location', 'SouthOutside');
    end
end

