function twoP_CSPcompare(area,expertise,depthRange)
% code to load imaging data from different animals and cell types to
% compare choice tuning of labeled vs non-labeled neurons.

%% some variables
rng(1);
allAnimals = {'CSP27' 'CSP29' 'CSP30'};

% expertise = 'expert';
% area = 'ALM';
% depthRange = [200 400];
dateRange = []; %use all if empty

% imgPath = 'G:\Google Drive\RawData\TwoP\richard_s2p_npy'; %path to imaging data
imgPath = '/Users/xiaonansun/Documents/twoP_data/';
% imgPath = 'X:\twoP'; %path to imaging data
bhvPath = 'G:\Google Drive\Behavior_Simon'; %path to behavioral data
sRate = 30.8987; %framerate of sutter 2p microscope
segIdx = [1 0.75 1.25 0.5 1];
segFrames = cumsum(floor(segIdx * sRate)); %max nr of frames per segment
groups = {'CSP'};
makePlot = false;
reload = false;

%% get data
clear allAuc allLeft allRight allRed allChoiceSelect
Cnt = zeros(1,3);
for iGroups = 1 : length(groups)
    
    cGroup = find(contains(allAnimals, groups{iGroups}))'; %mice in current group
    cMice = allAnimals(cGroup); %mice in current group
    Cnt = 0;
    
    for iAnimals = 1 : length(cMice)
        
        recNames = twoP_getRecordings(cMice{iAnimals}, expertise, depthRange, area, dateRange); % get recordings
        useMice(iAnimals) = ~isempty(recNames);
        
        for iRecs = 1 : length(recNames)
            Cnt = Cnt + 1;
            
            % get imaging and behavioral data
            baseFileName = [cMice{iAnimals} '_' recNames{iRecs}];
            disp(baseFileName);
            [Vc, redIdx, bhv] = twoP_checkImgBhvData_SM(cMice{iAnimals},recNames{iRecs}, imgPath, bhvPath, reload); %% Make sure imaging and behavior data are usable.
            [cAUC, cSelect, leftCells, rightCells] = twoP_checkChoiceSelectivity_SM(Vc, bhv, redIdx, segFrames, 70:100, makePlot);
            
            allAuc{iGroups, Cnt} = cAUC;
            allChoiceSelect{iGroups, Cnt} = cSelect;
            allLeft{iGroups, Cnt} = leftCells;
            allRight{iGroups, Cnt} = rightCells;
            allRed{iGroups, Cnt} = redIdx;
            drawnow;
            
        end
    end
end

[~,b] = cellfun(@size,allAuc);
fprintf('Mean number of cells per recording: %f, +- %f\n', nanmean(b), sem(b));

%% show PSTHs
figure('renderer','painters', 'name', sprintf('%s - depth: %d-%d, %s\n', area, depthRange(1), depthRange(2)));
clear cAUCs a b respFrac nrCells
groupLabels = groups;
groupLabels{end+1} = 'Unlabeled';
nrCols = 3;

for iGroups = 1 : length(groups)+1
    
    if iGroups == length(groups)+1 || ~isempty(cat(1, allLeft{iGroups,:}))
        
        if iGroups < length(groups)+1
            xx = 1 : size(allLeft,2);
            leftDat = cat(1, allLeft{iGroups,xx}); % concatednated single-cell neural activity averaged for all left-choice trials
            rightDat = cat(1, allRight{iGroups,xx}); % concatednated single-cell neural activity averaged for all right-choice trials
            cRed = logical(cat(1, allRed{iGroups,xx})); % concatenated indices of tdTomato-expressing neurons
            cSelect = cat(2, allChoiceSelect{iGroups,xx})'; % concatenated indices of choice-selective neurons, as defined as the top 5th percentile (2.5% on each end of the gaussian tail) of neurons based on the AUC analysis
            nrCells(iGroups) = sum(cRed);
            
            a = cat(2, allAuc{iGroups,xx}); % concatenated single-cell AUC values 
            cAUCs{iGroups} = a(cRed & cSelect);
            
            leftDat = leftDat(cRed & cSelect, :); % select tdTomato-expressing neurons that are left-choice-selective
            rightDat = rightDat(cRed & cSelect, :); % select tdTomato-expressing neurons that are right-choice-selective
            nanRej = isnan(mean(leftDat, 1)) | isnan(mean(rightDat, 1));
            leftDat(:, nanRej) = NaN;
            rightDat(:, nanRej) = NaN;
            
            a = [leftDat(:, ~nanRej), rightDat(:, ~nanRej)];
            a = zscore(a, [], 2);
            
%             a = a - nanmean(a(:, [1:15 sum(~nanRej) + (1:15)]),2);
%             baseMean = nanmean([leftDat(:, 1:10) rightDat(:, 1:10)], 2);
%             baseStd = nanstd([leftDat(:, 1:136) rightDat(:, 1:136)], [], 2);
%             a = bsxfun(@minus, a, baseMean);
%             a = bsxfun(@rdivide, a, baseStd);

            leftDat(:, ~nanRej) = a(:, 1 : size(a,2)/2);
            rightDat(:, ~nanRej) = a(:, size(a,2)/2 + 1 : end);
   
        else
            leftDat = cat(1, allLeft{xx});
            rightDat = cat(1, allRight{xx});
            cRed = logical(cat(1, allRed{xx}));
            cSelect = logical(cat(2, allChoiceSelect{xx}))';
            nrCells(iGroups) = sum(~cRed);

            leftDat = leftDat(~cRed & cSelect, :);
            rightDat = rightDat(~cRed & cSelect, :);
            nanRej = isnan(mean(leftDat, 1)) | isnan(mean(rightDat, 1));
            leftDat(:, nanRej) = NaN;
            rightDat(:, nanRej) = NaN;
            
            a = cat(2, leftDat(:, ~nanRej), rightDat(:, ~nanRej));
            a = zscore(a, [], 2);
            leftDat(:, ~nanRej) = a(:, 1 : size(a,2)/2);
            rightDat(:, ~nanRej) = a(:, size(a,2)/2 + 1 : end);
            
            cAUCs{iGroups} = cat(2, allAuc{xx});
            cAUCs{iGroups} = cAUCs{iGroups}(~cRed & cSelect);
        end
        
        % show total PSTH
        F = 1/sRate : 1/sRate : 136 / sRate;
        subplot(length(groups)+1,nrCols, ((iGroups-1)*nrCols) + 1);
        cLine(1) = stdshade(rightDat, 0.5, 'r', F); hold on;
        cLine(2) = stdshade(leftDat, 0.5, 'b', F); axis square
        title(groupLabels{iGroups}); xlabel('time(s)'); legend(cLine, {'ipsi', 'contra'}, 'location', 'northwest');
        ylabel('Event rates (SDUs)'); ylim([-1 2.5]); xlim([0 F(end)]);
        betterFigure;
        nvline(F([segFrames(1:2) segFrames(2)+30 segFrames(3)]+1), '--k')

        
        % show AUCs
        subplot(length(groups)+1,nrCols, ((iGroups-1)*nrCols) + 2);
        histogram(cAUCs{iGroups}, 0:0.04:1, 'DisplayStyle','stairs', 'Normalization', 'probability'); xlim([0 1]); axis square;
        title([groupLabels{iGroups} ' - AUC _C_h_o_i_c_e (selective only)']); ylabel('fraction of neurons');
        betterFigure;

%         % show Ipsi-contra distribution
        respFrac(1, iGroups) = sum(cAUCs{iGroups} < 0.5) / nrCells(iGroups); %ipsi response
        respFrac(2, iGroups) = sum(cAUCs{iGroups} > 0.5) / nrCells(iGroups); %contra response
        respFrac(3, iGroups) = 1 - sum(respFrac(1:2, iGroups)); %ipsi response
%         subplot(length(groups)+1,nrCols, ((iGroups-1)*nrCols) + 3);
%         pie(respFrac(:,iGroups), {'ipsi' 'contra' 'non-selective'})
%         title('Percent selective neurons')
%         betterFigure;

        % compare labeled vs unlabeled ispi/contra responses
        if iGroups == length(groups)+1
            
%             %unlabeled
%             leftDat = cat(1, allLeft{:});
%             rightDat = cat(1, allRight{:});
%             leftDat = leftDat(~cRed & cSelect, :);
%             rightDat = rightDat(~cRed & cSelect, :);
%             nanRej = isnan(mean(leftDat, 1)) | isnan(mean(rightDat, 1));
%             leftDat(:, nanRej) = NaN;
%             rightDat(:, nanRej) = NaN;
%             a = cat(2, leftDat(:, ~nanRej), rightDat(:, ~nanRej));
%             a = zscore(a, [], 2);
%             leftDat(:, ~nanRej) = a(:, 1 : size(a,2)/2);
%             rightDat(:, ~nanRej) = a(:, size(a,2)/2 + 1 : end);
%             
%             subplot(length(groups)+1,nrCols, 4);
%             line1(1) = stdshade(rightDat, 0.5, 'k', F); hold on; %ipsi response
%                        
%             %labeled
%             leftDat = cat(1, allLeft{:});
%             rightDat = cat(1, allRight{:});
%             leftDat = leftDat(cRed & cSelect, :);
%             rightDat = rightDat(cRed & cSelect, :);
%             nanRej = isnan(mean(leftDat, 1)) | isnan(mean(rightDat, 1));
%             leftDat(:, nanRej) = NaN;
%             rightDat(:, nanRej) = NaN;
%             a = cat(2, leftDat(:, ~nanRej), rightDat(:, ~nanRej));
%             a = zscore(a, [], 2);
%             leftDat(:, ~nanRej) = a(:, 1 : size(a,2)/2);
%             rightDat(:, ~nanRej) = a(:, size(a,2)/2 + 1 : end);
%             
%             subplot(length(groups)+1,nrCols, 4);
%             line1(2) = stdshade(rightDat, 0.5, 'g', F); hold on; %ipsi response
%             ylabel('Event rates (AUs)'); xlabel('time (s)');
%             title('IPSI responses'); legend(line1, 'All','CSP', 'location', 'northwest');
%             betterFigure;
%             
%             subplot(length(groups)+1,nrCols, ((iGroups-1)*nrCols) + 4);
%             line2(2) = stdshade(leftDat, 0.5, 'g', F); hold on; %contra response
%             ylabel('Event rates (AUs)'); xlabel('time (s)');
%             title('CONTRA responses'); legend(line2, 'All','CSP', 'location', 'northwest');
%             betterFigure;

            subplot(length(groups)+1,nrCols, 3);

            histogram(cAUCs{1}, 0:0.04:1, 'DisplayStyle','stairs', 'Normalization', 'probability', 'EdgeColor', 'k'); hold on;
            histogram(cAUCs{end}, 0:0.04:1, 'DisplayStyle','stairs', 'Normalization', 'probability', 'EdgeColor', 'g'); xlim([0 1]); axis square;

            xlim([0 1]); axis square;
            title('Distribution choice-selective neurons')
            betterFigure; axis square;
            

            subplot(length(groups)+1,nrCols, ((iGroups-1)*nrCols) + 3);
            respCells = round(respFrac .* nrCells); %total number of responsive neurons per group
            
            cCnt = 0.5; clear p
            for x = 1 : 3
                clear a b c
                cProb = respCells(x,end) / nrCells(end); %fraction of unlabeled neurons for current condition
                for y = 1 : length(groups)+1
                    [a(y),b(y),c(y)] = Behavior_wilsonError(respCells(x,y)*2,nrCells(y)*2); %error
                    
                    p(x,y) = myBinomTest(respCells(x,y),nrCells(y), cProb);
                    
                end
                bar(cCnt+1:cCnt+y,c); hold on;
                errorbar(cCnt+1:cCnt+y,c, a-c, c-b, '.k');
                
                cCnt = cCnt + y +2;
            end
            
%             [~, inGroupP(iGroups)] = binomCompare(respCells(1,iGroups), sum(respCells(1:2,iGroups)), ...
%                 respCells(2, iGroups), sum(respCells(1:2,iGroups)));
                        
            ax = gca;
            ax.XTick = (y+1)/2 : cCnt /3 : cCnt;
            ax.XTickLabel = {'ipsi' 'contra' 'non-selective'};
            title('Distribution choice-selective neurons')
            betterFigure; axis square;
            ax.XLim(2) = ax.XLim(2) + 1;
            
            for yy = 1 : iGroups
                    inGroupP(yy) = myBinomTest(respCells(1,yy), sum(respCells(1:2,yy)), respCells(2,yy)/sum(respCells(1:2,yy)));
            end
        end
    end
                
%             [~, inGroupP(iGroups)] = binomCompare(respCells(1,iGroups), sum(respCells(1:2,iGroups)), ...
%                 respCells(2, iGroups), sum(respCells(1:2,iGroups)));
end
fprintf('%s - depth: %d-%d, %s\n', area, depthRange(1), depthRange(2));
fprintf('P-choice, CSPvsAll: %d\n', ranksum(cAUCs{1}, cAUCs{length(groups)+1}));
fprintf('P-responseFraction, CSPvsAll: ipsi:%f, contra:%f, unresponsive:%f\n', p(1,1),p(2,1),p(3,1));
ipsiIdx = cAUCs{1} < 0.5;
fprintf('P-choice, CSP_IpsiVsContra: %f\n', ranksum(1-cAUCs{1}(ipsiIdx), cAUCs{1}(~ipsiIdx)));
fprintf('h for Fraction-choice, CSP_IpsiVsContra: %f\n', inGroupP(1));
fprintf('h for Fraction-choice, nonCSP_IpsiVsContra: %f\n', inGroupP(2));
fprintf('Labeled neurons:%d/%d (%f percent). %d recordings from %d mice.\n', sum(cRed), ...
    length(cRed), sum(cRed) / length(cRed),  Cnt, sum(useMice));

