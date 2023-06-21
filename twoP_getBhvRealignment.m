function newData = twoP_getBhvRealignment(data, cBhv, segFrames, opts)
% code to re-align data so each trial is aligned to different task episodes.
% Alignment is done to baseline, handle, stimulus and delay period.
% 'segFrames' defines how many frames per task episodes should be in the
% newData output. SegFrames should be cumulative, so say 'baseline should be
% from trial 1 to segFrames(1), handle data should be from segFrames(1) to
% segFrames(2) and so on...

%% align imaging data using Session data
rejCnt = 0;
newData = NaN(size(data,1), segFrames(5), size(data,3), 'single'); %new data to capture max duration of each segment
for iTrials = 1 : size(data,3)
    
    % get indices for current trial
    cIdx = 0;
    try
        stimOn = cBhv.RawEvents.Trial{iTrials}.Events.Wire3High; %time of stimulus onset - measured from soundcard
        handleOn = [reshape(cBhv.RawEvents.Trial{iTrials}.States.WaitForAnimal1',1,[]) ...
            reshape(cBhv.RawEvents.Trial{iTrials}.States.WaitForAnimal2',1,[]) ...
            reshape(cBhv.RawEvents.Trial{iTrials}.States.WaitForAnimal3',1,[])];

        clear cIdx
        cIdx(1) = handleOn(find(handleOn == cBhv.RawEvents.Trial{iTrials}.States.WaitForCam(1))-1); %find start of lever state that triggered stimulus onset
        cIdx(2) = stimOn;
        cIdx(3) = max(cat(2,cBhv.stimEvents{iTrials}{:})) + stimOn; %time of last stimulus event
        cIdx(4) = cBhv.RawEvents.Trial{iTrials}.States.MoveSpout(1);
        cIdx = floor((cIdx - stimOn + opts.preStim) * opts.frameRate); %convert to frames. This is the last frame of each segment.
        cIdx(end + 1) = size(data,2);
    end
    
    if cIdx(1) > 0 %in very rare cases there might be something wrong with handle time. don't use those trials.
        if segFrames(1) >= cIdx(1)
            newData(:, segFrames(1) - cIdx(1) + 1 : segFrames(1), iTrials) = data(:, 1 : cIdx(1), iTrials); % baseline
        elseif segFrames(1) < cIdx(1)
            newData(:, 1 : segFrames(1), iTrials) = data(:, cIdx(1) - segFrames(1) + 1 : cIdx(1), iTrials); % baseline
        end
           
        newData(:, segFrames(1) + 1 : segFrames(1) + (diff(cIdx(1:2))), iTrials) = data(:, cIdx(1) + 1 : cIdx(2), iTrials); %handle period
        newData(:, segFrames(2) + 1 : segFrames(2) + (diff(cIdx(2:3))), iTrials) = data(:, cIdx(2) + 1 : cIdx(3), iTrials); %stimulus period
        
        maxDiff = min([segFrames(3) + (diff(cIdx(3:4))) segFrames(5)]) - segFrames(3); %maximal possible delay duration
        newData(:, segFrames(3) + 1 : segFrames(3) + maxDiff, iTrials) = data(:, cIdx(3) + 1 : cIdx(3) + maxDiff, iTrials); %delay period
        
        if segFrames(4) + (diff(cIdx(4:5))) > segFrames(5)
            newData(:, segFrames(4) + 1 : segFrames(5), iTrials) = data(:, cIdx(4) + 1 : cIdx(4) + (segFrames(5) - segFrames(4)), iTrials); %response period
        else
            newData(:, segFrames(4) + 1 : segFrames(4) + (diff(cIdx(4:5))), iTrials) = data(:, cIdx(4) + 1 : cIdx(5), iTrials); %response period
        end
    else
        rejCnt = rejCnt + 1;
    end
end

if rejCnt > 0
    warning(['!!! Couldnt use ' num2str(rejCnt) ' trials because of broken handle initialization time or missing stimulus !!!'])
end