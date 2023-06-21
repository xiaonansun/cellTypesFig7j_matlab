function data = twoP_selectTwoPTrials(data,trials)
% Function to select a subset of trials/settings from selected 'trials' in a 
% larger array 'data' that has two-photon data. 'trials' should be a vector of 
% trial numbers that can be used.
% Usage: newData = twoP_selectTwoPTrials(data,trials)

%% get fieldnames
if isempty(data)
    bFields = {};
else
    bFields = fieldnames(data);
end

%% check if trials is logical index. If not create index that matches trials in data.
if ~islogical(trials)
    temp = false(1, length(data.trialNumbers));
    if any(trials > length(data.trialNumbers))
        warning('Trial index contains more trials as available in the dataset')
        trials(trials > length(data.trialNumbers)) = [];
    end
    temp(trials) = true;
    trials = temp;
end

%% cycle trough fields and carry over selected trials / sessions
for iFields = 1:size(bFields,1)
    if ~any(ismember(size(data.(bFields{iFields})), length(trials))) %if field does not contain single trials, it should contain session data instead
        if isstruct(data.(bFields{iFields})) %if field is a struct, check one layer deeper if it contains trial info
            tFields = fieldnames(data.(bFields{iFields}));
            if length(data.(bFields{iFields}).(tFields{1})) == length(trials)
                data.(bFields{iFields}).(tFields{1}) = data.(bFields{iFields}).(tFields{1})(trials);
            else
                data.(bFields{iFields}) =  data.(bFields{iFields}); %carry over complete field
            end
        else
            data.(bFields{iFields}) =  data.(bFields{iFields}); %carry over complete field
        end
    else
        if isvector(data.(bFields{iFields}))
            data.(bFields{iFields}) = data.(bFields{iFields})(trials); %carry over selected trials
        
        else %some highD matrix, find trial dimension, cut trials and reshape to match original matrix
            cIdx = find(ismember(size(data.(bFields{iFields})), length(trials))); %find trial dimension
            data.(bFields{iFields}) = arrayIndex(data.(bFields{iFields}), trials, cIdx); %get index from target dimension
                  
        end
    end 
end
