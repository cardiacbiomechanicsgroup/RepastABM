function [] = plotRepastABMStats(figureNum, varargin)
% This function plots the collagen and cell time course data from selected
% WoundABM model runs in the specified figure.

% INPUT:
%   figureNum: the integer figure number to plot into
%   varargin: structure containing the model parameters (structure),
%   statistics (structure), and collagen histograms (matrix).

% Author: Arlynn C Baker
% Created: 2019/06/30

numModels=nargin-1;

if numModels<1
    fprintf("plotRepastABMStats: No Models Specified");
else
    
    color=['r','g','b','k'];
    line={'-','--','-.',':'};
    
    figure(figureNum); hold on;
    for i=1:numModels
        
        % Determine plotting characteristics from model parameters
        mechs=varargin{i}.parameters.Mechanics;
        if isequal(mechs,'UniaxC')
            c=1;
        elseif isequal(mechs,'UniaxL')
            c=2;
        elseif isequal(mechs,'Biax')
            c=3;
        else
            c=4;
        end
        
        gridSize=varargin{i}.parameters.GridSize;
        if gridSize==10
            l=1;
        elseif gridSize==5
            l=2;
        elseif gridSize==2.5
            l=3;
        else
            l=4;
        end
        
        % Plot
        t=varargin{i}.statistics.Time;
        
        subplot(2,3,1); hold on;
        plot(t/24,varargin{i}.statistics.ColMVL,'Color',color(c),'LineStyle',line{l},'LineWidth',2);
        
        subplot(2,3,2); hold on;
        plot(t/24,varargin{i}.statistics.ColMVA,'Color',color(c),'LineStyle',line{l},'LineWidth',2);
        
        subplot(2,3,3); hold on;
        plot(t/24,varargin{i}.statistics.ColFRC,'Color',color(c),'LineStyle',line{l},'LineWidth',2);
        
        subplot(2,3,4); hold on;
        plot(t/24,varargin{i}.statistics.CellMVL,'Color',color(c),'LineStyle',line{l},'LineWidth',2);
        
        subplot(2,3,5); hold on;
        plot(t/24,varargin{i}.statistics.CellMVA,'Color',color(c),'LineStyle',line{l},'LineWidth',2);
        
        subplot(2,3,6); hold on;
        plot(t/24,varargin{i}.statistics.CellFRC,'Color',color(c),'LineStyle',line{l},'LineWidth',2);
    end
end
end