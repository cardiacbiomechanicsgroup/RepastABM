function [modelData] = readModelData(modelID)
% This function reads the statistics, histogram, and parameter output files
% of one or more WoundABM model runs.

% INPUT:
%   modelID: string containing the path information and model tag
%   (identifiers and time stamp) of the model. 
%   Ex: modelID='output\Biax_2.5_20190620_082011'

% OUTPUT: 
%   modelData: structure containing the model parameters (structure), 
%   statistics (structure), and collagen histograms (matrix).

% Author: Arlynn C Baker
% Created 2019/06/28
% Modified 2020/10/12

try
    modelData.statistics=readtable([modelID,'_WoundStat.csv']);
    rawHistograms=csvread([modelID,'_WoundColFiberAngDist.csv'],0,1);
    modelData.histograms=rawHistograms./sum((rawHistograms),2);
    modelData.parameters=table2struct(readtable([modelID,'_WoundParameters.csv']));
catch
    modelData.statistics=readtable([modelID,'_WoundStat_1.csv']);
    rawHistograms=csvread([modelID,'_WoundColFiberAngDist_1.csv'],0,1);
    modelData.histograms=rawHistograms./sum((rawHistograms),2);
    modelData.parameters=table2struct(readtable([modelID,'_WoundParameters_1.csv']));
end
try
    colMVL2D=csvread([modelID,'_colMVL2D.csv']);
    modelData.colMVL2D=colMVL2D(:,1:end-1);
    colMVA2D=csvread([modelID,'_colMVA2D.csv']);
    modelData.colMVA2D=colMVA2D(:,1:end-1);
    colFrac2D=csvread([modelID,'_colFrac2D.csv']);
    modelData.colFrac2D=colFrac2D(:,1:end-1);
catch
    try
    colMVL2D=csvread([modelID,'_colMVL2D_1.csv']);
    modelData.colMVL2D=colMVL2D(:,1:end-1);
    colMVA2D=csvread([modelID,'_colMVA2D_1.csv']);
    modelData.colMVA2D=colMVA2D(:,1:end-1);
    colFrac2D=csvread([modelID,'_colFrac2D_1.csv']);
    modelData.colFrac2D=colFrac2D(:,1:end-1);
    catch
        disp('no collagen heatmaps available');
    end
end
end