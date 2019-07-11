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

modelData.statistics=readtable([modelID,'_WoundStat.csv']);

rawHistograms=csvread([modelID,'_WoundColFiberAngDist.csv'],0,1);
modelData.histograms=rawHistograms./sum((rawHistograms),2);

modelData.parameters=table2struct(readtable([modelID,'_WoundParameters.csv']));
end