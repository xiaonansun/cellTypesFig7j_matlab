animal = 'CSP27';
bhvFileName = 'cBhv.mat';

clear bhvFilePathList

baseDir = '/Users/xiaonansun/Documents/data/twoP'; % change this directory as needed

% This generates file paths organized into a struct
bhvFilePathListStruct = dir(fullfile(baseDir,animal,'**',bhvFileName));

% This generates a list of behavior file paths organized as a cell
bhvFilePathList = cell(length(bhvFilePathListStruct),1);
for i = 1:length(bhvFilePathListStruct)
bhvFilePathList{i} = fullfile(bhvFilePathListStruct(i).folder,bhvFilePathListStruct(i).name);
end

%% Plots behavior performance of individual sessions

behavior_plotAllDiscCurves_RS(animal,bhvFilePathListStruct)

%%
bhv=behavior_LoadAllSessions_RS(animal,[],bhvFilePathListStruct);
behavior_PlotPerformanceMatrix(animal,bhv);
