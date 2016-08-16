function plot_electrode_statistic(stat, stat_name)    

n_electrodes = length(stat);

% plot median reliability across stimuli
n_electrodes_per_row = 50;
total_num_rows = ceil(n_electrodes / n_electrodes_per_row);
figure;
set(gcf,'Position',[0 0 1200 700]);
subplot(total_num_rows,1,1);
bar(1:n_electrodes,  stat, 'FaceColor', [0.5 0.5 0.5]);
set(gca, 'FontSize', 6);
xlim([0 n_electrodes+1]);
ylim([0 max(stat(:))*1.05]);
ylabel(stat_name);
for j = 1:ceil(n_electrodes / n_electrodes_per_row)
    subplot(total_num_rows+1,1,j+1);
    
    electrodes_to_plot = (1:n_electrodes_per_row) + n_electrodes_per_row*(j-1);
    xi = electrodes_to_plot<=n_electrodes;
    electrodes_to_plot = electrodes_to_plot(xi);
    
    bar(1:length(electrodes_to_plot),  stat(electrodes_to_plot), 'FaceColor', [0.5 0.5 0.5]);
    set(gca, 'FontSize', 6);
    xlim([0 n_electrodes_per_row+1]);
    ylim([0 max(stat(:))*1.05]);
    grid on;
    set(gca, 'XTick', 1:length(electrodes_to_plot), 'XTickLabel',electrodes_to_plot);
end

xlabel('Electrode');
