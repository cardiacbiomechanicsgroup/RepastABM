% This script plots the figures validating the Repast WoundABM against the 
% original MATLAB ABM (published J Physiol 2012) and experimental data 
% generated by Greg Fomovsky and Laura Caggiano.

% Author: Arlynn C Baker
% Created: 2020/09/17

clear; clc; close all;
addpath('C:\Users\acb4ft\Desktop\WoundABM\woundABMSim\output\Standard');
original_C=readModelData('UniaxC_10.0_20200917_161545');
original_L=readModelData('UniaxL_10.0_20200917_163634');
original_B=readModelData('Biax_10.0_20200917_160011');
original_models={original_C,original_L,original_B};

% Modified models
addpath('C:\Users\acb4ft\Desktop\WoundABM\woundABMSim\output\UniformNonCollagen_ScaleStructuralCue_NoRotation');
modified_B=readModelData('Biax_10.0_20201003_154812');
modified_C=readModelData('UniaxC_10.0_20201003_154812');
modified_L=readModelData('UniaxL_10.0_20201003_154812');
modified_models={modified_C,modified_L,modified_B};

% % NoFibrin_NoInitialColllagen_ScaleStructuralCue_NoRotation_MechanicalSweep
% addpath('C:\Users\acb4ft\Desktop\WoundABM\woundABMSim\output\NoNonCollagen_NoInitialCollagen_ScaleStructuralCue_NoRotation_MechanicsSweep');
% modified_C6=readModelData('UniaxC_10.0_20201007_130652_6');
% modified_L6=readModelData('UniaxL_10.0_20201007_130652_6');
% modified_B6=readModelData('Biax_10.0_20201007_130652_6');
% modified_models6={modified_C6,modified_L6,modified_B6};
% modified_C7=readModelData('UniaxC_10.0_20201007_131157_7');
% modified_L7=readModelData('UniaxL_10.0_20201007_131202_7');
% modified_B7=readModelData('Biax_10.0_20201007_131255_7');
% modified_models7={modified_C7,modified_L7,modified_B7};
% modified_C8=readModelData('UniaxC_10.0_20201007_131656_8');
% modified_L8=readModelData('UniaxL_10.0_20201007_131659_8');
% modified_B8=readModelData('Biax_10.0_20201007_131902_8');
% modified_models8={modified_C8,modified_L8,modified_B8};

figure(1); hold on;
styles={'r--','g--','b--'};
% styles2={'r:','g:','b:'};
% styles3={'r-.','g-.','b-.'};
% styles4={'r--','g--','b--'};
mechanics={'UniaxC','UniaxL','Biax'};
set(gcf, 'Position',  [0, 0, 1000, 300])

for i=1:3
    % Plotting
    subplot(1,3,i); hold on;
    t=modified_models{i}.statistics.Time;
    plot(t/24,modified_models{i}.statistics.ColMVL,styles{i},'LineWidth',1.5);
%     t=modified_models6{i}.statistics.Time;
%     plot(t/24,modified_models6{i}.statistics.ColMVL,styles2{i},'LineWidth',1.5);
%     plot(t/24,modified_models7{i}.statistics.ColMVL,styles3{i},'LineWidth',1.5);
%     plot(t/24,modified_models8{i}.statistics.ColMVL,styles4{i},'LineWidth',1.5);
    if strcmp(mechanics{i},'UniaxC')
        plot(t/24,original_models{i}.statistics.ColMVL,'r-','LineWidth',1.5);
        errorbar(21.538,0.433,(0.592-0.433),'o','Color','k','MarkerFaceColor','r','LineWidth',1.2);
    elseif strcmp(mechanics{i},'Biax')
        plot(t/24,original_models{i}.statistics.ColMVL,'b-','LineWidth',1.5);
        errorbar(20.551,0.103,(0.28-0.103),'o','Color','k','MarkerFaceColor','b','LineWidth',1.2);
    elseif strcmp(mechanics{i},'UniaxL')
        plot(t/24,original_models{i}.statistics.ColMVL,'g-','LineWidth',1.5);
        errorbar(7,0.52,(0.75-0.52),'o','Color','k','MarkerFaceColor','g','LineWidth',1.2);
        errorbar(14,0.39,(0.52-0.39),'o','Color','k','MarkerFaceColor','g','LineWidth',1.2);
        errorbar(21,0.38,(0.49-0.38),'o','Color','k','MarkerFaceColor','g','LineWidth',1.2);
        errorbar(42,0.42,(0.66-0.42),'o','Color','k','MarkerFaceColor','g','LineWidth',1.2);
    end
    
    % Formatting
    ylim([0 0.8]); yticks([0 0.2 0.4 0.6 0.8]);
    xticks([0 7 14 21 28 35 42]);
    xlim([0 42]);
    set(gca,'XTickLabel',{'0','7','14','21','28','35','42'},'FontName', 'Arial', 'FontSize', 12);
    set(gca, 'PlotBoxAspectRatio', [1.5 1 1], 'LineWidth', 1, 'FontName', 'Arial', 'FontSize', 12);
    if i==1
%         title('Circumferential Loading','FontName', 'Arial', 'FontSize', 12);
        ylabel('Collagen MVL', 'FontName', 'Arial', 'FontSize', 12);
    elseif i==2
%         title('Longitudinal Loading','FontName', 'Arial', 'FontSize', 12);
        xlabel('Time (d)', 'FontName', 'Arial', 'FontSize', 12);
    elseif i==3
%         title('Biaxial Loading','FontName', 'Arial', 'FontSize', 12);
        L(1) = plot(nan,nan,'-','Color',[.5,.5,.5],'LineWidth',1.5);
        L(2) = plot(nan(10,1),nan(10,1),':','Color',[.5,.5,.5],'LineWidth',1.5);
        L(3) = plot(nan,nan,'-o','Color','k','MarkerFaceColor',[.5,.5,.5],'LineWidth',1.5);
        legend(L,{'Standard','Modified','Experimental'},'Location','NorthEast','FontName','Arial','FontSize',12);
%         L(1) = plot(nan(10,1),nan(10,1),':','Color',[.5,.5,.5],'LineWidth',1.5);        
%         L(2) = plot(nan(10,1),nan(10,1),'-.','Color',[.5,.5,.5],'LineWidth',1.5);
%         L(3) = plot(nan(10,1),nan(10,1),'--','Color',[.5,.5,.5],'LineWidth',1.5);
%         legend(L,{'W_m=0.6','W_m=0.7','W_m=0.8'},'Location','NorthEast','FontName','Arial','FontSize',10);
    end
end