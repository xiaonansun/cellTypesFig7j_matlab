function [allAuc, cellSelect, leftCells, rightCells] = twoP_checkChoiceSelectivity_SM(Vc, bhv, redIdx, segFrames, timeIdx, makePlot)

recName = [bhv.Settings.SubjectName ', ' datestr(bhv.TrialStartTime(1), 'DD-MM-YYYY')];
plotCols = 4;
plotRows = 2;
% timeIdx = 70:110; %this is time during the trial from which to find AUC tuning
% timeIdx = 60:120; %this is time during the trial from which to find AUC tuning
nanCnt = sum(isnan(Vc(1,:,:)),3); %count number of NaNs across trials
nanThresh = 150; %dont use frames with more NaNs over trials as this
nrShuffles = 100; %nr of shuffles for bootstrapping

%% make a PSTH for left vs right choices
leftIdx = bhv.ResponseSide == 1 & bhv.Rewarded; % all left chosen trials
rightIdx = bhv.ResponseSide == 2 & bhv.Rewarded; % all left chosen trials

leftCells = nanmean(Vc(:,:,leftIdx),3);
rightCells = nanmean(Vc(:,:,rightIdx),3);
leftCells(:, nanCnt > nanThresh) = NaN;
rightCells(:, nanCnt > nanThresh) = NaN;

if makePlot
    figure('renderer', 'painters');
    
    % show PSTHs for left and right
    subplot(plotRows, plotCols, 1)
    leftMean = nanmean(leftCells(~redIdx, :),1);
    rightMean = nanmean(rightCells(~redIdx, :),1);
    
    a(1) = plot(leftMean, 'b'); hold on;
    a(2) = plot(rightMean, 'r'); axis square;
    betterFigure;
    title([recName ', non-redCells']);
    nvline(segFrames + 1,'--'); drawnow;
    legend(a, {'Left' 'Right'}, 'location', 'northwest');
    
    subplot(plotRows, plotCols, plotCols + 1)
    leftMean = nanmean(leftCells(redIdx, :),1);
    rightMean = nanmean(rightCells(redIdx, :),1);
    a(1) = plot(leftMean, 'b'); hold on;
    a(2) = plot(rightMean, 'r'); axis square;
    betterFigure;
    nvline(segFrames + 1,'--'); drawnow;
    title([recName ', redCells']); xlabel('time(frames)');
    legend(a, {'Left' 'Right'}, 'location', 'northwest');
end


%% check for tuned neurons in the defined time range
rightIdx = bhv.ResponseSide == 2; % all left chosen trials
rightIdx = rightIdx(bhv.Rewarded);

X = squeeze(nanmean(Vc(:, timeIdx, bhv.Rewarded), 2)); %this is neurons by trials

%compute AUC for all neurons
allAuc = colAUC(X', rightIdx, 'abs', false);
testAuc = zeros(nrShuffles, size(X,1), 'single');
for x = 1 : nrShuffles
    testAuc(x, :) = colAUC(X', rightIdx(randperm(length(rightIdx))), 'abs', false);
end

% identify choice selective cells
cellSelect = allAuc > prctile(testAuc,97.5,1) | allAuc < prctile(testAuc,2.5,1);
    
if makePlot
    subplot(plotRows, plotCols, 2)
    histogram(allAuc(~redIdx),50); xlim([0 1]); axis square
    betterFigure;
    title('Choice AUC, non-redCells');
    
    subplot(plotRows, plotCols, plotCols + 2)
    histogram(allAuc(redIdx),50); xlim([0 1]); axis square
    betterFigure;
    title('Choice AUC, redCells');
    
end

%% single cell examples
if makePlot
    
    leftIdx = bhv.ResponseSide == 1; % all left chosen trials
    rightIdx = bhv.ResponseSide == 2; % all right chosen trials

    % ipsi-selective non-labeled example
    subplot(plotRows, plotCols, 3)
    cIdx = find(~redIdx);
    [~, b] = min(allAuc(~redIdx));
    cData = squeeze(Vc(cIdx(b), :, :));
    cData(nanCnt > nanThresh, :) = NaN;
    
    stdshade(cData(:, leftIdx)', 0.5, 'b'); hold on;
    stdshade(cData(:, rightIdx)', 0.5, 'r');
    nvline(segFrames + 1,'--'); axis square; drawnow;
    title('non-labeled, ipsi example');
    
    % contra-selective non-labeled example
    subplot(plotRows, plotCols, 4)
    cIdx = find(~redIdx);
    [~, b] = max(allAuc(~redIdx));
    cData = squeeze(Vc(cIdx(b), :, :));
    cData(nanCnt > nanThresh, :) = NaN;
    
    stdshade(cData(:, leftIdx)', 0.5, 'b'); hold on;
    stdshade(cData(:, rightIdx)', 0.5, 'r');
    nvline(segFrames + 1,'--'); axis square; drawnow;
    title('non-labeled, contra example');
    
    
    % ipsi-selective non-labeled example
    subplot(plotRows, plotCols, plotCols + 3)
    cIdx = find(redIdx);
    [~, b] = min(allAuc(redIdx));
    cData = squeeze(Vc(cIdx(b), :, :));
    cData(nanCnt > nanThresh, :) = NaN;
    
    stdshade(cData(:, leftIdx)', 0.5, 'b'); hold on;
    stdshade(cData(:, rightIdx)', 0.5, 'r');
    nvline(segFrames + 1,'--'); axis square; drawnow;
    title('red ipsi example');
    
    % contra-selective non-labeled example
    subplot(plotRows, plotCols, plotCols + 4)
    cIdx = find(redIdx);
    [~, b] = max(allAuc(redIdx));
    cData = squeeze(Vc(cIdx(b), :, :));
    cData(nanCnt > nanThresh, :) = NaN;
    
    stdshade(cData(:, leftIdx)', 0.5, 'b'); hold on;
    stdshade(cData(:, rightIdx)', 0.5, 'r');
    nvline(segFrames + 1,'--'); axis square; drawnow;
    title('red contra example');
end

%% identify selective neurons based on ROC at all times
% if makePlot
%
%     nrFrames = size(Vc,2);
%     leftIdx = bhv.ResponseSide == 1; % all left chosen trials
%     rightIdx = bhv.ResponseSide == 2; % all left chosen trials
%     stepSize = 1;
%
%     allAuc = NaN(size(Vc,1), floor(nrFrames/stepSize));
%     Cnt = 0;
%     for x = 1 : stepSize : (nrFrames - stepSize)+1
%
%         Cnt = Cnt + 1;
%
%         % isolate the right frames for the current timepoint
%         X = squeeze(nanmean(Vc(:, x : x + stepSize-1, :), 2)); %this is neurons by trials for the current time point
%         useIdx = ~isnan(bhv.ResponseSide) & ~isnan(X(1,:)); %use only trials with a response
%
%         if sum(useIdx) > 50 %needs a minimum of usable trials
%
%             X = X(:, useIdx); %only use valid trials for data
%             Y = leftIdx(useIdx); %only use valid trials for choice variable
%
%             %compute AUC for all neurons
%             allAuc(:,Cnt) = colAUC(X', Y, 'abs', false);
%
%         end
%     end
%
%
%     %show auc matrix
%     subplot(plotRows, plotCols, 2)
%     imagesc(allAuc(~redIdx,:)); axis image;
%     try colormap(colormap_BlueWhiteRed(256)); end
%     caxis([0.25 0.75])
%     title('auc values for current recording'); colorbar;
%
%     subplot(plotRows, plotCols, plotCols + 2)
%     imagesc(allAuc(redIdx,:)); axis image;
%     try colormap(colormap_BlueWhiteRed(256)); end
%     caxis([0.25 0.75])
%     title('auc values for current recording'); colorbar;
% end
