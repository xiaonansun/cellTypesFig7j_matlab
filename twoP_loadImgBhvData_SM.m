function [npy,data,SessionData,bhvFilePath,suite2pDir]=twoP_loadImgBhvData_SM(varargin)

animal = varargin{1};
session = varargin{2};
imagingRootDir = varargin{3};
bhvRootDir = varargin{4};

if nargin > 5
    doSmooth = varargin{5};
    filterWindow = varargin{6};
else
    doSmooth = true;
    filterWindow = 10;
end

analysisFileName = [animal '_' session];

imagingSubDir = 'suite2p\plane0';

disp(['Animal: ' animal '; Session: ' session]);

% --- Load the google sheets document "2photon acquisition record" --- %
docid = '16MKB18byS2S7ATSopf2NJpVPFZnH-BJdh4NtvBY80_k'; %Irene's copy of google sheet
expTable=GetGoogleSpreadsheet(docid); % this function (GetGoogleSpreadsheet.m) needs to be downloaded
bhvColIdx=find(contains(expTable(1,:),'Behavior file name'));
iFolderColIdx=find(contains(expTable(1,:),'Folder'));
rejFramesColIdx=find(contains(expTable(1,:),'discard frames'));

try
    bhvRowIdx = find(ismember(expTable(:,iFolderColIdx),session));
    bhvFName = expTable{bhvRowIdx(ismember(expTable(bhvRowIdx),animal)),bhvColIdx};
    rejFrames = expTable{bhvRowIdx(ismember(expTable(bhvRowIdx),animal)),rejFramesColIdx};
    try
        rejFrames = str2num(rejFrames);
    catch
        rejFrames = [];
    end
catch ME
    disp([ME.message]);
    disp('Cannot load session. Please check the session name input.');
    analysis.error.behaviorTable = ME.message;
end


fprintf('Loading 2P imaging data...');

% organize suite2p data into npy files (ROI and spiking)
tic
npy = twoP_importSuite2pData(animal,session, imagingRootDir);
disp(['The .npy file was loaded in ' num2str(toc) ' seconds.'])
npy.ops.rejFrames = rejFrames;

% IMPORTANT!!! smoothes inferred spikes with a gaussian filter
if doSmooth==true
    npy.spks = smoothCol(npy.spks,2,filterWindow,'gauss');
end

tic
data = twoP_alignDetectionTask_SM(npy.ops, npy, npy.iscell, npy.redcell, npy.bin_MScan_filepath); % align suite2p data to sensory stimulus
disp(['Trial alignment to stimulus onset was completed in ' num2str(toc) ' seconds.'])

if exist('ME','var')
    disp(['Error detected: ' ME.identifier '. ' ME.message])
else
    disp('2P DATA LOADED!');
end

% Load behavior data
fprintf('Loading behavior data...');
[SessionData,bhvFilePath] = twoP_loadBehaviorSession_IL(animal,session,bhvFName,bhvRootDir); 
fprintf('DONE!\n');

suite2pDir = fullfile(imagingRootDir,animal,'imaging',session,imagingSubDir); data.suite2pDir = suite2pDir;
disp(['Directory of suite2p output: ' suite2pDir]);
disp(['Path of behavior data: ' bhvFilePath]);
disp(['Number of trial codes received by MScan: ' num2str(max(data.trialNumbers))]);
disp(['Number of trials included in neural data matrix (MScan analog data not rejected): ' num2str(length(data.trialNumbers))]);
disp(['Number of trials recorded by Bpod: ' num2str(length(SessionData.Rewarded))]);
if abs(data.maxTrialCnt - length(SessionData.Rewarded)) >= 10; disp('Behavior and imaging differs by more than 10 trials: do you have the correct behavior file?'); end

end