function [outliers, multiscale_outliers] ...
    = envelope_outliers(envelopes, env_sr, figure_directory, varargin)

% Detects outliers in signal envelopes. The difference between the median and
% the 84th percential of the envelope distribution for each electrode is
% measured. For a Gaussian distribution this interval is equal to a standard
% deviation. Envelopes are considered outliers if they fall a certain number of
% 'standard deviations' above the median. 
% 
% 2016-08-15 - Created, Sam NH

I.scales = [0 0.25 1 4];
I.thresholds = [5 2.5 2 1.5];
I.plot_mosaic = true;
I.plot_individ_electrodes = false;
I = parse_optInputs_keyvalue(varargin, I);

% total number of channels
n_channels = size(envelopes,2);

% make scale/threshold parameters electrode specific if not already
if isvector(I.scales);
    I.scales = repmat(I.scales(:), 1, n_channels);
else
    assert(size(I.scales) == n_channels);
end
if isvector(I.thresholds);
    I.thresholds = repmat(I.thresholds(:), 1, n_channels);
else
    assert(size(I.thresholds) == n_channels);
end
assert(all(size(I.thresholds)==size(I.scales)));
n_scales = size(I.scales,1);

% normalize by 90% of the distribution
s = quantile(envelopes,0.9);
Z = envelopes ./ repmat(s, size(envelopes,1),1);
clear sd;

% time vector used for plotting
t = (0:size(envelopes,1)-1)/env_sr;

if I.plot_individ_electrodes
    individual_electrode_directory = [figure_directory '/individual-electrodes'];
    if ~exist(individual_electrode_directory, 'dir');
        mkdir(individual_electrode_directory);
    end
end

%% loop through channels

if I.plot_mosaic || I.plot_individ_electrodes
    figh = figure;
end
fig_pdfs = cell(1,n_channels);
outliers = false(size(envelopes));
multiscale_outliers = false([size(envelopes), n_scales]);
for q = 1:n_channels
    
    if I.plot_individ_electrodes
        clf(figh);
        set(figh, 'Position', [100 100 1200 n_scales*150]);
    end
    
    for i = 1:n_scales
        if I.scales(i,q) == 0
            z = Z(:,q);
        else
            s = round(I.scales(i)*env_sr);
            h = normpdf(-3*s:3*s,0,s)';
            z = myconv(min(Z(:,q),10), h, 'causal', false, 'norm', true);
        end
        
        multiscale_outliers(:,q,i) = z > I.thresholds(i,q);
        
        if I.plot_individ_electrodes
            sbpt = subplot(n_scales,1,i);
            set(sbpt,'Position', [0.1, 1-(1/n_scales)*i+0.06, 0.8, 0.15]);
            plot(t, z, 'r-', 'LineWidth', 0.1); hold on;
            z(multiscale_outliers(:,q,i)) = NaN;
            plot(t, z, 'k-', 'LineWidth', 0.1);
            plot([0 t(end)], I.thresholds(i,q)*[1 1], 'r--', 'LineWidth', 0.1);
            xlim([0 t(end)]);
            ylim([0 I.thresholds(i,q)*3]);
        end
    end
    
    % save
    if I.plot_individ_electrodes
        % save as matlab figure
        saveas(gcf, [individual_electrode_directory '/electrode' num2str(q) '-allscales.fig']);
        
        % save as pdf
        set(gcf, 'PaperSize', [8, n_scales*2]);
        set(gcf, 'PaperPosition', [0.25, 0.25, 7.5, n_scales*2-0.5]);
        print([individual_electrode_directory '/electrode' num2str(q) '-allscales.pdf'],'-dpdf');
    end
    
    % collapse across all outliers
    outliers(:,q) = any(multiscale_outliers(:,q,:), 3);
    
    % plot with outliers from al scales
    if I.plot_individ_electrodes
        clf(figh);
        set(figh, 'Position', [100 100 1200 150]);
        z = Z(:,q);
        plot(t, z, 'r-', 'LineWidth', 0.1); hold on;
        z(outliers(:,q)) = NaN;
        plot(t, z, 'k-', 'LineWidth', 0.1);
        xlim([0 t(end)]);
        ylim([0 10]);
        title(sprintf('elec %d', q));
        
        % save as matlab figure
        saveas(gcf, [individual_electrode_directory '/electrode' num2str(q) '.fig']);
        
        % save as pdf
        set(gcf, 'PaperSize', [8, 2]);
        set(gcf, 'PaperPosition', [0.25, 0.25, 7.5, 1.5]);
        print([individual_electrode_directory '/electrode' num2str(q) '.pdf'],'-dpdf');
        
        % expand along x-dimension and save as pdf
        box off;
        sbpt = subplot(1,1,1);
        set(sbpt,'Position', [0 0 1 1]);
        set(gca, 'XTick',[],'YTick',[]);
        fig_pdfs{q} = [individual_electrode_directory '/electrode' num2str(q) '-expanded.pdf'];
        set(gcf, 'PaperSize', [40, 1]);
        set(gcf, 'PaperPosition', [1, 0.25, 38, 0.5]);
        print(fig_pdfs{q},'-dpdf');
    end        
end

% pdf with sets of electrodes for ease of viewing
if I.plot_individ_electrodes
    for i = 1:ceil(n_channels/20);
        xi = (1:20) + (i-1)*20;
        xi = intersect(xi, 1:size(envelopes,2));
        append_pdfs([figure_directory '/mosaic-electrodes' num2str(xi(1)) '-' num2str(xi(end)) '.pdf'], fig_pdfs{xi});
    end
end

if I.plot_mosaic
    
    mean_outliers = 100*mean(any(multiscale_outliers, 3), 1);
    plot_electrode_statistic(mean_outliers, '% Outliers');
    
    % save as pdf
    fig_dims = get(gcf, 'Position');
    set(gcf, 'PaperSize', fig_dims(3:4)/130);
    set(gcf, 'PaperPosition', [0.25, 0.25, (fig_dims(3:4)/130-0.25)]);
    print([figure_directory '/mean-number-of-outliers.pdf'],'-dpdf');
    saveas(gcf, [figure_directory '/mean-number-of-outliers.fig']);
    
end


% 
% q = 186;
% 
% outliers = false(size(env_z(:,q)));
% for i = 1:length(scales)
%     if scales(i) == 0
%         env_smoothed = env_z(:,q);
%     else
%         N = round(scales(i)*3);
%         h = normpdf(-N:N,0,scales(i))';
%         env_smoothed = myconv(env_z(:,q), h, 'causal', false, 'norm', true);
%     end
%     outliers(env_smoothed>thresholds(i)) = 1;
%     figure;
%     plot(env_smoothed)
% end
% 
% 
% % plot(env_smoothed); ylim([0 3]);
% 
% figure;
% z = env_z(:,q);
% z(outliers) = NaN;
% set(gcf, 'Position', [0 0 1400 300]);
% t = (0:size(env_z,1)-1)/env_sr;
% plot(t, env_z(:,q),'r'); hold on;
% plot(t,z,'b'); 
% xL = xlim;
% plot(xL, thresholds(1) * [1 1], 'r--');
% ylim([0 10]);
% 
% %%
% 
% 
% % plot
% plot_electrode_statistic(mean(outliers)*100, '% Outliers');
% box off;
% set(gcf, 'PaperSize', [12 6]);
% set(gcf, 'PaperPosition', [0.25 0.25 11.5 5.5]);
% print([figure_fname '.pdf'],'-dpdf');
% print([figure_fname '.png'],'-dpng', '-r100');
% close all;
% 
% % 
% figure;
% set(gcf, 'Position', [0 0 1400 300]);
% plot(env_z(:,80)); hold on;
% xL = xlim;
% plot(xL, outlier_threshold * [1 1], 'r--');
% ylim([0 10]);
% 
% %%
% 
% 
% %%
% 
% x = ones(10,1);
% h = ones(4,1);
% y = conv(x,h,'same');
% plot(y)
% 
% %%
% 
% x_pad = [zeros(size(h,1),1); x]
% h_pad = [h; zeros(size(x_pad,1)-size(h,1),1)]
% 
% y = ifft2(fft2(x_pad) .* fft2(h_pad));
% y = y(size(h,1)+1:end,:)
% 
% 
% 
% %%
% 
% 
% x = ones(size(envelopes,1),1);
% h = ones(200,1)/200;
% y = conv(x,h,'same');
% 
% global_outliers = conv(double(outliers(:,186)),h,'same');
% 
% 
% x = resample(env_z(:,186), 1, 100)
% s = diff(quantile(x,[0.5 0.9]));
% x = x ./ repmat(s, size(x,1),1);
% plot(x);
% ylim([0 10]);
% 
% N = 100;
% x = outliers(:,80);
% n_outliers = sum(x(1:N*2+1));
% global_outliers = zeros(size(x));
% for i = N+1:length(x)-N;
%     n_outliers = n_outliers + x(i-1+N) - x(i-1);
%     global_outliers(i) = n_outliers;
% end