function [] = plotRepastABMCollagenHeatMaps(figureNum, type, varargin)
% This function plots the 2D collagen heatmaps and calculates the
% heterogeneity metrics for the Repast WoundABM.

% INPUT:
%   figureNum: the integer figure number to plot into
%   varargin: structure containing the model information

% Author: Arlynn C Baker
% Created: 2020/10/12

numModels=(nargin-2);

mapMVA = [0 1 0; 0.05 0.95 0; 0.1 0.9 0; 0.15 0.85 0; 0.2 0.8 0;
    0.25 0.75 0; 0.3 0.7 0; 0.35 0.65 0; 0.4 0.6 0; 0.45 0.55 0; 0.5 0.5 0;
    0.55 0.45 0; 0.6 0.4 0; 0.65 0.35 0; 0.7 0.3 0; 0.75 0.25 0; 0.8 0.2 0;
    0.85 0.15 0; 0.9 0.1 0; 0.95 0.05 0; 1 0 0; 0.95 0.05 0; 0.9 0.1 0;
    0.85 0.15 0; 0.8 0.2 0; 0.75 0.25 0; 0.7 0.3 0; 0.65 0.35 0; 0.6 0.4 0;
    0.55 0.45 0; 0.5 0.5 0; 0.45 0.55 0; 0.4 0.6 0; 0.35 0.65 0; 0.3 0.7 0;
    0.25 0.75 0; 0.2 0.8 0; 0.15 0.85 0; 0.1 0.9 0; 0.05 0.95 0; 0 1 0];

if numModels==0
    fprintf("plotRepastABMCollagenHeatMaps: No Models Specified");
else
    
    figure(figureNum); hold on;
    tiledlayout(2,3)
    for i=1:numModels
        
        % Plot heatmaps
        nexttile;
        if isequal(type,'MVL')
            imagesc(varargin{i}.colMVL2D);
            colormap(spring)
            title(['collisionGuidance=',num2str(varargin{i}.parameters.CollisionGuidance)]);
            
        elseif isequal(type,'MVA')
            imagesc(varargin{i}.colMVA2D);
            colormap(mapMVA)
            title(['collisionGuidance=',num2str(varargin{i}.parameters.CollisionGuidance)]);
            
        elseif isequal(type,'Frac')
            imagesc(varargin{i}.colFrac2D);
            colormap(spring)
            title(['collisionGuidance=',num2str(varargin{i}.parameters.CollisionGuidance)]);
        end
    end
    set(gcf, 'Position',  [0, 0, 1000, 600])
end
end