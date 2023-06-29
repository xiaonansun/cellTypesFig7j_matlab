animal = 'CSP27';
bhvFileName = 'cBhv.mat';

clear bhvFilePathList

baseDir = '/Users/xiaonansun/Documents/data/twoP';

bhvFilePathListStruct = dir(fullfile(baseDir,animal,'**',bhvFileName));
bhvFilePathList = cell(length(bhvFilePathListStruct),1);
for i = 1:length(bhvFilePathListStruct)
bhvFilePathList{i} = fullfile(bhvFilePathListStruct(i).folder,bhvFilePathListStruct(i).name);
end

%%

behavior_plotAllDiscCurves_RS(animal,bhvFilePathListStruct)