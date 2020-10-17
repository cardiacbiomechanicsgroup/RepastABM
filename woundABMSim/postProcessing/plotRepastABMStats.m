function [] = plotRepastABMStats(figureNum, styles, varargin)
% This function plots the collagen and cell time course data from selected
% WoundABM model runs in the specified figure with the specified styles.

% INPUT:
%   figureNum: the integer figure number to plot into
%   styles: array containing line styles
%   varargin: structure containing the model parameters (structure),
%   statistics (structure), and collagen histograms (matrix).

% Author: Arlynn C Baker
% Created: 2019/06/30
% Modified: 2020/02/17

numModels=(nargin-2);

if numModels==0
    fprintf("plotRepastABMStats: No Models Specified");
else

    figure(figureNum); hold on;
    for i=1:numModels
        
        % Plot
        t=varargin{i}.statistics.Time;
        
        subplot(2,3,1); hold on;
        plot(t/24,varargin{i}.statistics.ColMVL,styles{i},'LineWidth',2);
        
        subplot(2,3,2); hold on;
        plot(t/24,varargin{i}.statistics.ColMVA,styles{i},'LineWidth',2);
        
        subplot(2,3,3); hold on;
        plot(t/24,varargin{i}.statistics.ColFRC,styles{i},'LineWidth',2);
        
        subplot(2,3,4); hold on;
        plot(t/24,varargin{i}.statistics.CellMVL,styles{i},'LineWidth',2);
        
        subplot(2,3,5); hold on;
        plot(t/24,varargin{i}.statistics.CellMVA,styles{i},'LineWidth',2);
        
        subplot(2,3,6); hold on;
        plot(t/24,varargin{i}.statistics.CellFRC,styles{i},'LineWidth',2);
    end
end
end