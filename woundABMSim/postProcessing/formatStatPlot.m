function [] = formatStatPlot(figureNum, legendEntries)
% This function applies the figure formatting for the ABM cell and collagen
% time course statistics.

% INPUT:
%     figureNum: the integer figure number to plot into
%     legendEntries: an array of strings describing the models plotted in
%     the figure

% Author: Arlynn C Baker
% Created: 2019/06/21

figure(figureNum); hold on;

subplot(2,3,1); hold on;
ylabel('Collagen MVL', 'FontName', 'Arial', 'FontSize', 14);
ylim([0 0.8]); yticks([0 0.2 0.4 0.6 0.8]);
setTimeAxisLabels();

subplot(2,3,2); hold on;
ylabel('Collagen MVA', 'FontName', 'Arial', 'FontSize', 14);
ylim([-90 90]); yticks([-90 -60 -30 0 30 60 90]);
setTimeAxisLabels();

subplot(2,3,3); hold on;
ylabel('Collagen Fraction', 'FontName', 'Arial', 'FontSize', 14);
legend(legendEntries, 'Location', 'NorthWest');
ylim([0 0.5]); yticks([0 0.1 0.2 0.3 0.4 0.5]);
setTimeAxisLabels();

subplot(2,3,4); hold on;
ylabel('Cell MVL', 'FontName', 'Arial', 'FontSize', 14);
ylim([0 0.8]); yticks([0 0.2 0.4 0.6 0.8]);
setTimeAxisLabels();

subplot(2,3,5); hold on;
ylabel('Cell MVA', 'FontName', 'Arial', 'FontSize', 14);
ylim([-90 90]); yticks([-90 -60 -30 0 30 60 90]);
setTimeAxisLabels();

subplot(2,3,6); hold on;
ylabel('Cell Fraction', 'FontName', 'Arial', 'FontSize', 14);
ylim([0 1]); yticks([0 0.2 0.4 0.6 0.8 1.0]);
setTimeAxisLabels();
end

function [] = setTimeAxisLabels()
xlabel('Time (d)', 'FontName', 'Arial', 'FontSize', 14)
xticks([0 7 14 21 28 35 42])
xlim([0 42]);
set(gca,'XTickLabel',{'0','7','14','21','28','35','42'},'FontName', 'Arial', 'FontSize', 14);
set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'LineWidth', 1, 'FontName', 'Arial', 'FontSize', 14);
end