function [] = plotMatlabABMStats(figureNum, mechanics)
% This function plots the collagen and cell time course data from the
% MATLAB ABM (J Pysiol 2012).

% INPUTS:
%     figureNum: the integer figure number to plot into
%     mechanics: the model mechanics to be plotted ('UniaxC', 'UniaxL', 'Baix')

% Author: Arlynn C Baker
% Created: 2019/06/30
% Modified: 2020/02/17

% Load JPhysiology2012 ABM data
load('ABM_20130314_ResultsSummary.mat','T','CollTimeMVL','CollTimeMVA','CollTimeFRC','CellTimeMVL','CellTimeMVA','CellTimeFRC');

% Interpret mechanics
if isequal(mechanics,'Biax')
    type=1;
    style=[0,0,1];
elseif isequal(mechanics,'UniaxC')
    type=2;
    style=[1,0,0];
elseif isequal(mechanics,'UniaxL')
    type=3;
    style=[0,1,0];
end

figure(figureNum); hold on;

subplot(2,3,1); hold on;
plot(T, CollTimeMVL(:,type), 'Color', style, 'LineWidth',1.5);

subplot(2,3,2); hold on;
plot(T, CollTimeMVA(:,type), 'Color', style, 'LineWidth',1.5);

subplot(2,3,3); hold on;
plot(T, CollTimeFRC(:,type), 'Color', style, 'LineWidth',1.5);

subplot(2,3,4); hold on;
plot(T, CellTimeMVL(:,type), 'Color', style, 'LineWidth',1.5);

subplot(2,3,5); hold on;
plot(T, CellTimeMVA(:,type), 'Color', style, 'LineWidth',1.5);

subplot(2,3,6); hold on;
plot(T, CellTimeFRC(:,type), 'Color', style, 'LineWidth',1.5);
end