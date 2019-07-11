function [] = formatHistPlot(figureNum, legendEntries)
% This function applies the figure formatting for the ABM collagen
% histograms.

% INPUT:
%     figureNum: the integer figure number to plot into
%     legendEntries: an array of strings describing the models plotted in
%     the figure

% Author: Arlynn C Baker
% Created: 2019/06/30

figure(figureNum); hold on;
sgtitle('WoundABM: time-course histogram');
for j=1:1:6
    subplot(2,3,j); hold on;
    title(strcat('Week ', num2str(j)));
    xlabel('Collagen MVA (Degree)', 'FontName', 'Arial', 'FontSize', 12);
    ylabel('Collagen Fiber Fraction', 'FontName', 'Arial', 'FontSize', 12);
    xlim([-90,90]); xticks([-90 -60 -30 0 30 60 90]);
    ylim([0 0.1]); yticks([0 0.02 0.04 0.06 0.08 0.1]);
    set(gca,'XTickLabel',{'-90','-60','-30','0','30','60','90'},'FontName', 'Arial', 'FontSize', 12);
    set(gca,'PlotBoxAspectRatio', [1.5 1 1], 'LineWidth', 1, 'FontName','Arial', 'FontSize', 14);
    if j==3
        legend(legendEntries);
    end
end
set(gcf,'pos',[600 400 700 400]);
end