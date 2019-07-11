function [] = plotRepastABMHist(figureNum, varargin)
% This function plots the collagen histograms at each week in the selected 
% WoundABM model runs in the specified figure.

% INPUT:
%   figureNum: the integer figure number to plot into
%   varargin: structure containing the model parameters (structure),
%   statistics (structure), and collagen histograms (matrix).

% Author: Arlynn C Baker
% Created: 2019/06/30

numModels=nargin-1;

if numModels<1
    fprintf("plotRepastABMHist: No Models Specified");
else
    
    color=['r','g','b','k'];
    line={'-','--','-.',':'};
    
    angleBin=-87.5:5:87.5;
    
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
        for j=1:1:6
            subplot(2,3,j); hold on;
            plot(angleBin,varargin{i}.histograms(7*j,1:end),'Color',color(c),'LineStyle',line{l},'LineWidth',2);
        end
    end
end
end