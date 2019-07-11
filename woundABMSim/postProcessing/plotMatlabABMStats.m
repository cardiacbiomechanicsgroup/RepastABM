function [] = plotMatlabABMStats(figureNum, mechanics)
% This function plots the collagen and cell time course data from the
% MATLAB ABM (J Pysiol 2012), and Greg Fomovsky's experimental data onto
% the specified figure.

% INPUTS:
%     figureNum: the integer figure number to plot into
%     mechanics: the model mechanics to be plotted ('UniaxC', 'UniaxL', 'Baix')

% Author: Arlynn C Baker
% Created: 2019/06/30

% Load JPhysiology2012 ABM data
load('ABM_20130314_ResultsSummary.mat','T','CollTimeMVL','CollTimeMVA','CollTimeFRC','CellTimeMVL','CellTimeMVA','CellTimeFRC');

% Load Greg's data
load('FomovskyData.mat','ColExpMean','ColExpSD','ExpTime');

% Interpret mechanics
if isequal(mechanics,'Biax')
    type=1;
elseif isequal(mechanics,'UniaxC')
    type=2;
elseif isequal(mechanics,'UniaxL')
    type=3;
end

figure(figureNum); hold on;

subplot(2,3,1); hold on;
plot(T, CollTimeMVL(:,type), '-k', 'LineWidth', 2);
if strcmp(type,'UniaxC')
    errorbar(21.538,0.433,(0.592-0.433),'o','Color',[.5,.5,.5],'MarkerFaceColor',[.5,.5,.5],'LineWidth',1);
elseif strcmp(type,'Biax')
    errorbar(20.551,0.103,(0.28-0.103),'o','Color',[.5,.5,.5],'MarkerFaceColor',[.5,.5,.5],'LineWidth',1);
end

subplot(2,3,2); hold on;
plot(T, CollTimeMVA(:,type), '-k', 'LineWidth', 2);

subplot(2,3,3); hold on;
plot(T, CollTimeFRC(:,type), '-k', 'LineWidth', 2);
errorbar(ExpTime,ColExpMean,ColExpSD,'o','Color',[.5,.5,.5],'MarkerFaceColor',[.5,.5,.5],'LineWidth', 1);

subplot(2,3,4); hold on;
plot(T, CellTimeMVL(:,type), '-k', 'LineWidth', 2);

subplot(2,3,5); hold on;
plot(T, CellTimeMVA(:,type), '-k', 'LineWidth', 2);

subplot(2,3,6); hold on;
plot(T, CellTimeFRC(:,type), '-k', 'LineWidth', 2);
end