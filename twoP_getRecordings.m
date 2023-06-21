function recNames = twoP_getRecordings(animal, expertise, depthRange, area, dateRange)
% function to select specific recordings in an animal of interest.
% Based on animalID, expertise, cortical depth, area, and training time.

if ~exist('dateRange','var')
    dateRange = [];
end

%% --- Load the google sheets document "2photon acquisition record" --- %
docid = '16MKB18byS2S7ATSopf2NJpVPFZnH-BJdh4NtvBY80_k'; %google sheet with recordings

expTable=GetGoogleSpreadsheet(docid); % this function (GetGoogleSpreadsheet.m) needs to be downloaded

animalColIdx = contains(expTable(1,:),'Animal');
expertiseColIdx = contains(expTable(1,:),'Expertise');
depthColIdx = contains(expTable(1,:),'Depth');
locationColIdx = contains(expTable(1,:),'Location');
dateColIdx = contains(expTable(1,:),'ExperimentDate');
folderColIdx = contains(expTable(1,:),'Folder');


%% make selection
animalRowIdx = find(strcmpi(expTable(:,animalColIdx),animal));

% select correct entries and adjust animalRowIdx
expertiseSelect = strcmpi(expTable(animalRowIdx,expertiseColIdx),expertise);

cDepths = cellfun(@(x) str2num(x), expTable(animalRowIdx,depthColIdx));
depthSelect = cDepths >= depthRange(1) & cDepths <= depthRange(2);

areaSelect = strcmpi(expTable(animalRowIdx,locationColIdx),area);

if ~isempty(dateRange) && iscell(dateRange)
    dateRange = cellfun(@(x) datenum(x), dateRange);
    cDates = cellfun(@(x) datenum(x), expTable(animalRowIdx,dateColIdx));
    dateSelect = cDates >= dateRange(1) & cDates <= dateRange(2);
else
    dateSelect = true(length(animalRowIdx), 1);
end

% correct seelction
cSelect = animalRowIdx(expertiseSelect & depthSelect & areaSelect & dateSelect);
recNames = expTable(cSelect,folderColIdx);

