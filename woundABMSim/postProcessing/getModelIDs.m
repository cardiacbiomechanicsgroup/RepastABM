function [modelIDs] = getModelIDs()
% This function retrieves the unique file ID's of the WoundABM model
% outputs from the current directory.

% OUTPUT: fileIDs is an array of the unique output file ID's in the current
%   directory. The unique file ID does not contian the output type 
%   identifier (_WoudParameters.csv, _WoundColFibAngDist.csv, or 
%   _WoundStat.csv).

% Author: Arlynn C. Baker
% Created: 2019/07/24

folderInfo=dir('*_WoundParameters*.csv');
modelIDs=cell(1);
for i=1:length(folderInfo)
    fileName=folderInfo(i).name;
    idIndex=strfind(fileName,'_WoundParameters');
    modelIDs{i}=fileName(1:idIndex(end)-1);
end
end