function plot_all_envelopes(exp, subjid, r, MAT_file_envelope)

% general-purpose ecog analysis code
global root_directory;

% directory for this project
project_directory = [root_directory '/' exp];

load(MAT_file_envelope, 'envelopes', 'env_sr', 'band_in_Hz');

% MAT file to save results to
bp_freq_range_string = ...
    [num2str(band_in_Hz(1)) '-' num2str(band_in_Hz(2)) 'Hz'];

% directory to save plots to
figure_directory = [project_directory '/figures/envelopes/' ...
    subjid '/r' num2str(r) '/bpfilt_' bp_freq_range_string];
if ~exist(figure_directory, 'dir');
    mkdir(figure_directory);
end

sd = quantile(envelopes,0.9);
envelopes = envelopes ./ repmat(sd, size(envelopes,1),1);

figh = figure;
fig_pdfs = cell(1,size(envelopes,2));
for i = 1:size(envelopes,2)
   
    t = (1:size(envelopes,1))/env_sr;
    clf(figh);
    set(figh, 'Position', [100 100 1200 200]);
    sbpt = subplot(1,1,1);
    set(sbpt,'Position', [0 0 1 1]);
    plot(t, envelopes(:,i), 'k-', 'LineWidth', 0.1);
    xlim([0 t(end)]); ylim([0 4]);
    set(gca, 'YTick', []);
    title(sprintf('elec %d',i));
    
    % save as matlab figure
    saveas(gcf, [figure_directory '/electrode' num2str(i) '.fig']);    
    
    % save as pdf
    fig_pdfs{i} = [figure_directory '/electrode' num2str(i) '.pdf'];
    set(gcf, 'PaperSize', [40, 1]);
    set(gcf, 'PaperPosition', [1, 0.25, 38, 0.5]);
    print(fig_pdfs{i},'-dpdf');
    
end

for i = 1:ceil(size(envelopes,2)/20);
    xi = (1:20) + (i-1)*20;
    xi = intersect(xi, 1:size(envelopes,2));
    append_pdfs([figure_directory '/mosaic-electrodes' num2str(xi(1)) '-' num2str(xi(end)) '.pdf'], fig_pdfs{xi});
end

