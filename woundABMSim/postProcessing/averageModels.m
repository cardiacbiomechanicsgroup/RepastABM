function [averageModelData] = averageModels(varargin)
% This function averages multiple WoundABM models. Before averaging, the 
% function will check the model parameters to ensure that they are 
% compatible. Note that all model statistics must contain the same 
% statistics, time points, and parameters. 

% INPUT: 
%   varargin: structure containing the model parameters (structure), 
%   statistics (structure), and collagen histograms (matrix).

% OUTPUT:
%   averageModelData: structure containing the average model parameters
%   (structure), statistics (structure), and collagen histograms (matrix).

% Author: Arlynn C baker
% Created: 2019/06/29

numModels=nargin;
if numModels==0
    fprintf("averageModels Error: No Models Were Specified \n");
elseif numModels==1
    averageModelData=varargin{1};
else
    test=0;
    for i=1:numModels
    paramNames{i}=fieldnames(varargin{i}.parameters);
    end
    for i=2:numModels
        test=test+not(isequal(paramNames{i-1},paramNames{i}));
    end
    if test==0
        statNames=varargin{1}.statistics.Properties.VariableNames;
        numStats=length(statNames);
        numEntries=height(varargin{1}.statistics);
        statistics=table;
        statSum=cell(1:numStats);
        histSum=0;
        for j=1:numStats
            statSum{j}=zeros(numEntries,1);
            for i=1:numModels
                statSum{j}=statSum{j}+getfield(varargin{i}.statistics,statNames{j});
            end
            statistics=setfield(statistics,statNames{j},statSum{j}/numModels);
        end
        for i=1:numModels
            histSum=histSum+varargin{i}.histograms;
        end
        averageModelData.statistics=statistics;
        averageModelData.histograms=histSum/numModels;
        averageModelData.parameters=varargin{1}.parameters;
    end
end
end