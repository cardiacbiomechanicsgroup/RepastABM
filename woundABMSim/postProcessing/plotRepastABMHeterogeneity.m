function [] = plotRepastABMHeterogeneity(figureNum, varargin)
% This function plots the fiber-fiber dot products vs fiber-fiber distances

% INPUT:
%   figureNum: the integer figure number to plot into
%   varargin: fiberInfoStruct's containing the aggregated model information

% Author: Arlynn C Baker
% Created: 2020/15/12

numModels=(nargin-1);

if numModels==0
    fprintf("plotRepastABMCollagenHeatMaps: No Models Specified");
else
    
    figure(figureNum); hold on;
    if numModels==1
        ylabel('Fiber-Fiber Alignment', 'FontName', 'Arial', 'FontSize', 12);
    end
    t=tiledlayout(2,3);
    for i=1:numModels
        nexttile; hold on;
        
        [AveFiberDotDistances,AveFiberDotMVL]=heterogeneityAnalysis(varargin{i});
        
        plot(AveFiberDotDistances(:,end),AveFiberDotMVL(:,end),'LineWidth', 1.2)
        title(['CollisionGuidance=',num2str(varargin{i}.CollisionGuidance)])
        ylim([0.5 1]); xlim([0 220])
        set(gca,'FontName', 'Arial', 'FontSize', 12);
        set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'LineWidth', 1, 'FontName', 'Arial', 'FontSize', 12);
        
    end
    xlabel(t,'Fiber-Fiber Distance (\mum)')
    ylabel(t,'Fiber-Fiber Alignment')
    set(gcf, 'Position',  [0, 0, 1000, 500])
    
    figure(figureNum+1); hold on;
    for i=1:numModels
        [AveFiberDotDistances,AveFiberDotMVL]=heterogeneityAnalysis(varargin{i});
        if i==1
            plot(AveFiberDotDistances(:,end),AveFiberDotMVL(:,end),'Color',[0,0,0],'LineWidth', 1.2)     
        else
            plot(AveFiberDotDistances(:,end),AveFiberDotMVL(:,end),'Color',[1,(i-2)/(numModels),(i-2)/(numModels)],'LineWidth', 1.2)     
        end
    end
end
set(gcf, 'Position',  [0, 0, 500, 300])
legend('Standard','CG=0.05','CG=0.1','CG=0.2','CG=0.4','CG=0.8');
ylim([0.5 1]); xlim([0 220])
set(gca,'FontName', 'Arial', 'FontSize', 12);
set(gca, 'LineWidth', 1, 'FontName', 'Arial', 'FontSize', 12);
xlabel('Fiber-Fiber Distance (\mum)', 'FontName', 'Arial', 'FontSize', 12)
ylabel('Fiber-Fiber Alignment', 'FontName', 'Arial', 'FontSize', 12);
end