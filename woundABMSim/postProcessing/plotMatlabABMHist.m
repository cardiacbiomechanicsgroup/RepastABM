function [] = plotMatlabABMHist(figureNum, mechanics)
% This function plots the collagen histograms at each week from the MATLAB 
% ABM (J Pysiol 2012), and Greg Fomovsky's experimental data onto
% the specified figure.

% INPUTS:
%     figureNum: the integer figure number to plot into
%     mechanics: the model mechanics to be plotted ('UniaxC', 'UniaxL', 'Baix')

% Author: Arlynn C Baker
% Created: 2019/06/30

% Load JPhysiology2012 ABM data
load('ABM_20130314_ResultsSummary.mat','FiberDirBins','CollHist3WK');

% Interpret mechanics
if isequal(mechanics,'Biax')
    type=1;
elseif isequal(mechanics,'UniaxC')
    type=2;
elseif isequal(mechanics,'UniaxL')
    type=3;
end

figure(figureNum); hold on;
subplot(2,3,3); hold on;
plot(FiberDirBins, CollHist3WK(:,type), '-k', 'LineWidth', 2);
end