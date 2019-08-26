function [] = plotBatchData( modelIDs, stats, batchParams )
% This function plots the time course data for a specific output statistics
% from the WoundABM model runs specified with the titles

% cd \\bme-jwh-nas02.bme.virginia.edu\Home\ColRotationSweep_batch
% addpath('C:\Users\abake\Desktop\WoundABM\woundABMSim\postProcessing');
% modelIDs=getModelIDs;
% stats={'ColMVL','ColMVA','ColFRC','CellMVL','CellMVA','CellFRC'};
% batchParams={'Wp','Wm','Wc','Ws'};
% plotBatchData(modelIDs,stats,batchParams);

% INPUT:
%   modelIDs: array containing the model IDs
%   params: array containing the output parameters (strings) to plot
%   batchSpecs: array containing the input parameters (strings) to title
%   the subplots by

% Author: Arlynn C Baker
% Created: 2019/07/24

color=['r','g','b','k'];
line={'-','--','-.',':'};

figureNum=1;
numModels=length(modelIDs);
subDim=ceil(sqrt(numModels));

% Iterate over parameters
for i=1:length(stats)
    
    % Set figure number
    figure(figureNum);
    if ishandle(figureNum)
        figureNum=figureNum+1;
    end
    
    % Iterate over models
    for j=1:numModels
        
        % Read in model data
        modelData=readModelData(modelIDs{j});
        
        % Determine plotting characteristics from model parameters
        mechs=modelData.parameters.Mechanics;
        if isequal(mechs,'UniaxC')
            c=1;
        elseif isequal(mechs,'UniaxL')
            c=2;
        elseif isequal(mechs,'Biax')
            c=3;
        else
            c=4;
        end
        gridSize=modelData.parameters.GridSize;
        if gridSize==10
            l=1;
        elseif gridSize==5
            l=2;
        elseif gridSize==2.5
            l=3;
        else
            l=4;
        end
        
        % Assemble plot title by iterating over batch parameters
        for k=1:length(batchParams)
            paramVal=modelData.parameters.(batchParams{k});
            if k==1
                params=[batchParams{k},'=',num2str(round(paramVal,2))];
            else 
                params=[params,', ',batchParams{k},'=',num2str(round(paramVal,2))];
            end
        end
        
        % Plot
        subplot(subDim,subDim,j);
        t=modelData.statistics.Time;
        plot(t/24,modelData.statistics.(stats{i}),'Color',color(c),...
            'LineStyle',line{l},'LineWidth',2);
        title(params);
    end
end
end