function data = twoP_combineStimAlignedData_SM(data)
% 2021-10-04 This script combines stimulus-aligned 2P data primary to deal
% with MScan computer crashes. When the MScan PC crashes, imaging is
% terminated while the behavior continues, hence one or more trials will be
% skipped in the 2P data. This results in missing trial codes from the
% analog input. The script twoP_alignDetectionTask.m has been modified to
% combine multiple 2P sub-sessions into a single data struct. This data
% struct separates individual sub-sessions into cells. This script
% concatenates data.trialNumbers, data.stimFrameOrig, data.neural, data.DS,
% and data.analog.

if iscell(data.neural)
    fprintf('Merging data: %s - %s\n', data.animal, data.session)
    data.analogFreq = data.analogFreq{1};
    data.trialNumbers = horzcat(data.trialNumbers{1},data.trialNumbers{2});
    data.trialStimFrame = data.trialStimFrame{1};
    data.stimFramesOrig = horzcat(data.stimFramesOrig{1},data.stimFramesOrig{2});
    data.neural = cat(3,data.neural{1},data.neural{2});
    data.neuralTimes = data.neuralTimes{1};
    data.neuralTimesMs = data.neuralTimesMs{1};
    data.DS = cat(3,data.DS{1},data.DS{2});
    data.trialStimSample = data.trialStimSample{1};
    data.stimSamplesOrig = horzcat(data.stimSamplesOrig{1},data.stimSamplesOrig{2});
    data.analog = cat(3,data.analog{1},data.analog{2});
    data.analogTimes = data.analogTimes{1};
end