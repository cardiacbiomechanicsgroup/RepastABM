function [] = plotExperimentalData(figureNum, mechanics)
% This function plots experimental collagen and cell time course onto the
% specified figure.

% INPUTS:
%     figureNum: the integer figure number to plot into
%     mechanics: the model mechanics to be plotted ('UniaxC', 'UniaxL', 'Baix')

% Author: Arlynn C Baker
% Created: 2020/02/17
% Modified: 2020/10/15

figure(figureNum); hold on;

% Collagen MVL
subplot(2,3,1); hold on;
if strcmp(mechanics,'UniaxC') % Greg's cryoinfarct data 2012
    errorbar(21.538,0.433,(0.592-0.433),'o','Color','k','MarkerFaceColor','r','LineWidth',1.2);
elseif strcmp(mechanics,'Biax') % Greg's cryoinfarct data 2012
    errorbar(20.551,0.103,(0.28-0.103),'o','Color','k','MarkerFaceColor','b','LineWidth',1.2);
elseif strcmp(mechanics,'UniaxL') % Laura's patch data 2018
    errorbar(7,0.52,(0.75-0.52),'o','Color','k','MarkerFaceColor','g','LineWidth',1.2);
    errorbar(14,0.39,(0.52-0.39),'o','Color','k','MarkerFaceColor','g','LineWidth',1.2);
    errorbar(21,0.38,(0.49-0.38),'o','Color','k','MarkerFaceColor','g','LineWidth',1.2);
    errorbar(42,0.42,(0.66-0.42),'o','Color','k','MarkerFaceColor','g','LineWidth',1.2);
end

% Collagen Fraction
subplot(2,3,3); hold on;
if strcmp(mechanics,'UniaxC') % Greg's cryoinfarct data 2012
    errorbar(21,0.349,0.07,'o','Color','k','MarkerFaceColor','r','LineWidth',1.2);
elseif strcmp(mechanics,'Biax') % Greg's cryoinfarct data 2012
    errorbar(21,0.349,0.067,'o','Color','k','MarkerFaceColor','b','LineWidth',1.2);
elseif strcmp(mechanics,'UniaxL') % Laura's patch data 2018
    errorbar(7,0.087,(0.14-0.087),'o','Color','k','MarkerFaceColor','g','LineWidth',1.2);
    errorbar(14,0.25,(0.35-0.25),'o','Color','k','MarkerFaceColor','g','LineWidth',1.2);
    errorbar(21,0.32,(0.50-0.32),'o','Color','k','MarkerFaceColor','g','LineWidth',1.2);
    errorbar(42,0.48,(0.60-0.48),'o','Color','k','MarkerFaceColor','g','LineWidth',1.2);
end

% % Greg's infarct data 2010 (used in original ABM paper)
% load('Fomovsky2010Data.mat','ColExpMean','ColExpSD','ExpTime');
% % Collagen Fraction
% subplot(2,3,3); hold on;
% if strcmp(mechanics,'UniaxC')
%     errorbar(ExpTime,ColExpMean,ColExpSD,'o','Color','k','MarkerFaceColor','r','LineWidth',1.2);
% elseif strcmp(mechanics,'Biax')
%     errorbar(ExpTime,ColExpMean,ColExpSD,'o','Color','k','MarkerFaceColor','b','LineWidth',1.2);
% elseif strcmp(mechanics,'UniaxL') 
%     errorbar(ExpTime,ColExpMean,ColExpSD,'o','Color','k','MarkerFaceColor','g','LineWidth',1.2);
% end
end