function [Vc, redIdx, cBhv, suite2p_output_dir, VcNeuropil]=twoP_checkImgBhvData_SM(varargin)

animal = varargin{1};
session = varargin{2};
imgPath = varargin{3};
bhvPath = varargin{4};
reload = varargin{5};

imagingSubDir = 'suite2p\plane0';

disp(['Animal: ' animal '; Session: ' session]);

fprintf('Loading 2P imaging data...\n');

% organize suite2p data into npy files (ROI and spiking)
tic

% suite2p_output_dir = [imgPath filesep animal filesep 'imaging' filesep session filesep 'suite2p\plane0'];
suite2p_output_dir = [imgPath filesep animal filesep 'SpatialDisc' filesep session];

if ~reload
    try
        
        load([suite2p_output_dir filesep 'Vc.mat'], 'Vc', 'redIdx');
        load([suite2p_output_dir filesep 'VcNeuropil.mat'], 'VcNeuropil', 'redIdx');
        load([suite2p_output_dir filesep 'cBhv.mat'], 'cBhv');
        if ~exist('redIdx', 'var')
            error;
        end
       
    catch
        reload = true;
    end
end

if reload
    
    [~,data,SessionData,~,~] = twoP_loadImgBhvData_SM(animal,session, imgPath, bhvPath, true, 10);
    data = twoP_combineStimAlignedData_SM(data); %check if data needs to be merged
    
    [rejIdx, ~] = twoP_checkTrialDrift(data, 0.4, 0.6, [], false); %check for drift over the session
    data.neural = data.neural(~rejIdx, :, :); %dont use drifty neurons
    data.dFOF = data.dFOF(~rejIdx, :, :); %dont use drifty neurons
    data.FNeu = data.FNeu(~rejIdx, :, :); %dont use drifty neurons
    data.idx_redcell = data.idx_redcell(~rejIdx);
    
    trialIdx = data.trialNumbers(data.trialNumbers <= length(SessionData.Rewarded));
    bhv = selectBehaviorTrials(SessionData, trialIdx); %subselect trials that are in the twoP data structure
    data = twoP_selectTwoPTrials(data, 1:length(trialIdx)); %make sure that twoP structure does not have more trials as the behavior
    
    %% aligne data to trial episodes
    opts.preStim = data.trialStimFrame*data.msPerFrame/1000; % Duration of the data (in seconds) before the stimulus occurs
    opts.frameRate = 1000/data.msPerFrame; % Frame rate of imaging
    segIdx = [1 0.75 1.25 0.5 1];
    segFrames = cumsum(floor(segIdx * opts.frameRate)); %max nr of frames per segment
    
    redIdx = data.idx_redcell;
    Vc = twoP_getBhvRealignment(data.neural, bhv, segFrames, opts); %aligned to different trial episodes
    VcNeuropil = twoP_getBhvRealignment(data.FNeu, bhv, segFrames, opts); %aligned to different trial episodes
    cBhv = bhv;
    
    save([suite2p_output_dir filesep 'data.mat'], 'data');
    save([suite2p_output_dir filesep 'Vc.mat'], 'Vc', 'redIdx');
    save([suite2p_output_dir filesep 'VcNeuropil.mat'], 'VcNeuropil', 'redIdx');
    save([suite2p_output_dir filesep 'cBhv.mat'], 'cBhv');

end