% This script plots the 2D collagen heatmaps and calculates the
% heterogeneity metrics for the Repast WoundABM.

% Author: Arlynn C Baker
% Created: 2020/10/12
clear; clc; close all;

addpath('C:\Users\acb4ft\Desktop\WoundABM\woundABMSim\output\CollisionGuidance_6Week');

% Models with no initial structure, no rotation, a scaled structural cue,
% Wm=0.6, and pastMigration is inhibited at 6 weeks
% biax_0_v1=readModelData('Biax_10.0_20201016_003918_0');
% biax_05_v1=readModelData('Biax_10.0_20201016_000959_05');
% biax_1_v1=readModelData('Biax_10.0_20201015_235515_1');
% biax_2_v1=readModelData('Biax_10.0_20201016_001830_2');
% biax_4_v1=readModelData('Biax_10.0_20201016_002513_4');
% biax_8_v1=readModelData('Biax_10.0_20201016_003231_8');

% Models with no initial structure, no rotation, a scaled structural cue,
% Wm=0.6, and pastMigration is inhibited at 2 weeks
biax_0_v1=readModelData('Biax_10.0_20201016_220022_0');
biax_05_v1=readModelData('Biax_10.0_20201016_220022_05');
biax_1_v1=readModelData('Biax_10.0_20201016_220022_1');
biax_2_v1=readModelData('Biax_10.0_20201016_220151_2');
biax_4_v1=readModelData('Biax_10.0_20201016_220154_4');
biax_8_v1=readModelData('Biax_10.0_20201016_220151_8');

% Plot heat maps
plotRepastABMCollagenHeatMaps(1, 'MVA', biax_0_v1, biax_05_v1, biax_1_v1, biax_2_v1, biax_4_v1, biax_8_v1);

% Consolidate multiple model runs
[biaxCombined_0]=packageRepastABMHeterogeneityData(biax_0_v1);%, biax_0_v2, biax_0_v3, biax_0_v4);
[biaxCombined_05]=packageRepastABMHeterogeneityData(biax_05_v1);%, biax_05_v2);%, biax_05_v3, biax_05_v4);
[biaxCombined_1]=packageRepastABMHeterogeneityData(biax_1_v1);%, biax_1_v2, biax_1_v3, biax_1_v4);
[biaxCombined_2]=packageRepastABMHeterogeneityData(biax_2_v1);%, biax_2_v2, biax_2_v3, biax_2_v4);
[biaxCombined_4]=packageRepastABMHeterogeneityData(biax_4_v1);%, biax_4_v2, biax_4_v3, biax_4_v4);
[biaxCombined_8]=packageRepastABMHeterogeneityData(biax_8_v1);%, biax_8_v2, biax_8_v3, biax_8_v4);

% Plot heterogeneity
plotRepastABMHeterogeneity(2,biaxCombined_0,biaxCombined_05,biaxCombined_1,biaxCombined_2,biaxCombined_4,biaxCombined_8);