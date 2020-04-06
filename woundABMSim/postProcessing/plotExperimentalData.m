function [] = plotExperimentalData(figureNum, mechanics)
% This function plots experimental collagen and cell time course onto the
% specified figure.

% INPUTS:
%     figureNum: the integer figure number to plot into
%     mechanics: the model mechanics to be plotted ('UniaxC', 'UniaxL', 'Baix')

% Author: Arlynn C Baker
% Created: 2020/02/17

% Load Greg's data
load('FomovskyData.mat','ColExpMean','ColExpSD','ExpTime');

figure(figureNum); hold on;

subplot(2,3,1); hold on;
if strcmp(mechanics,'UniaxC')
    errorbar(21.538,0.433,(0.592-0.433),'o','Color',[.5,.5,.5],'MarkerFaceColor',[.5,.5,.5],'LineWidth',1);
elseif strcmp(mechanics,'Biax')
    errorbar(20.551,0.103,(0.28-0.103),'o','Color',[.5,.5,.5],'MarkerFaceColor',[.5,.5,.5],'LineWidth',1);
elseif strcmp(mechanics,'UniaxL')   % Estimated
    fprintf('Estimated Longitudinal DataSet in plotExperimentalData');
    errorbar(7,0.52,(0.75-0.52),'o','Color',[.5,.5,.5],'MarkerFaceColor',[.5,.5,.5],'LineWidth',1);
    errorbar(14,0.39,(0.52-0.39),'o','Color',[.5,.5,.5],'MarkerFaceColor',[.5,.5,.5],'LineWidth',1);
    errorbar(21,0.38,(0.49-0.38),'o','Color',[.5,.5,.5],'MarkerFaceColor',[.5,.5,.5],'LineWidth',1);
    errorbar(42,0.42,(0.66-0.42),'o','Color',[.5,.5,.5],'MarkerFaceColor',[.5,.5,.5],'LineWidth',1);
end

subplot(2,3,3); hold on;
errorbar(ExpTime,ColExpMean,ColExpSD,'o','Color',[.5,.5,.5],'MarkerFaceColor',[.5,.5,.5],'LineWidth', 1);
end