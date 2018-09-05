function allData = prepTemporalData(allData, varargin)
%% Convert snapshot data into temporal trajectories

[   minCount                                            , ...
    noMissVar                                           , ...
    staggered                                           , ...
    resample                                            , ...
    samplingRate            ] = process_options(varargin, ...
    'minCount'              , 2                         , ...
    'noMissVar'             , []                        , ...
    'staggered'             , false                     , ...
    'resample'              , false                     , ...
    'samplingRate'          , 6                         );

data            = allData.data;
isSorted        = isDataSorted(data);

assert(isSorted, 'Data table is not sorted!');

% Figure out which patients have missing data amongst the subset
% 'noMissVar'
if ~isempty(noMissVar)
    phiMask             = false(size(allData.phiNames));
    for p=1:length(noMissVar)
        phiMask         = phiMask | strcmpi(allData.phiNames, noMissVar{p});
    end
    subPhi              = allData.phi(:, phiMask);
    mask                = ~any(isnan(subPhi), 2);
    
    [patients, ~, uRid]     = unique(data.RID);
    counts                  = histc(uRid, 1:length(patients));
    patientMask             = mat2cell(mask, counts);
    patientMask             = cellfun(@(seq)sum(seq) >= minCount, ...
        patientMask);
    mask                    = patientMask(uRid);
    
    data                    = data(mask, :);
    allData.phi             = allData.phi(mask, :);
    allData.profile         = allData.profile(mask, :);
end

% Figure out which patients have at least minCount occurences
rid                     = data.RID;
[patients, ~, ridIdx]   = unique(rid);
ridCounts               = histc(ridIdx, 1:length(patients));
ridCounts               = ridCounts(ridIdx);
mask                    = ridCounts >= minCount;

data                    = data(mask, :);
allData.phi             = allData.phi(mask, :);
allData.profile         = allData.profile(mask, :);
allData.data            = data;

[isSorted, ~, order]    = isDataSorted(data);
if ~isSorted
    warning('Data not sorted!');
    allData.phi             = allData.phi(order, :);
    allData.profile         = allData.profile(order, :);
end

if resample
    allData                 = padData(allData, samplingRate);
end

data                        = allData.data;
rid                         = data.RID;
[~, lastIdx, ~]             = unique(rid, 'last');

compMat         = bsxfun(@eq, colvec(rid), rowvec(rid));

trajectoryIndices       = cell(length(rid), 1);

for ii=1:length(rid)
    trajectoryIndices{ii}   = colvec(find(compMat(ii, 1:ii)));
end

if ~staggered
    trajectoryIndices       = trajectoryIndices(lastIdx);
end

allData.trajIdx         = trajectoryIndices;
allData.counts          = cellfun(@length, trajectoryIndices);

end

function allData = padData(allData, res)

data            = allData.data;
rid             = data.RID;
patients        = unique(rid);

linInterPhi     = {'CONVTIME_L', 'CONVTIME_U', 'CMCI_L', 'CMCI_U', ...
    'FOLLOWUP', 'AGE'};
holdInterPhi    = {'RID', 'CAT_GENDER_MALE', 'CAT_GENDER_FEMALE', ...
    'CAT_APGEN1_2_0', 'CAT_APGEN1_3_0', 'CAT_APGEN1_4_0', ...
    'CAT_APGEN2_2_0', 'CAT_APGEN2_3_0', 'CAT_APGEN2_4_0', ...
    'PTEDUCAT', 'PTWORK', 'PTETHCAT', 'PTRACCAT', 'ANARTERR'};
dataColNames    = data.Properties.VariableNames;

linDataColIdx   = zeros(length(linInterPhi), 1);
holdDataColIdx  = zeros(length(holdInterPhi), 1);

% the variables that will be linearly interpolated
for ii=1:length(linInterPhi)
    linDataColIdx(ii)      = find(strcmpi(dataColNames, linInterPhi{ii}));
end
% the variables that we apply zero-order hold to
for ii=1:length(holdInterPhi)
    holdDataColIdx(ii)      = find(strcmpi(dataColNames, holdInterPhi{ii}));
end
% the DX variable is treated specially
diagnosisIdx    = find(strcmpi(dataColNames, 'DX'));

% the numeric visit codes
numViscode      = getNumericViscode(data.VISCODE2);
paddedData      = cell(length(patients), 1);
paddedPhi       = cell(length(patients), 1);

for p=1:length(patients)
    patientMask         = rid == patients(p);
    visTimes            = round(numViscode(patientMask)/res)*res;
    samplingTimes       = diff(visTimes);
    
    patData             = table2cell(data(patientMask, :));
    patPhi              = allData.phi(patientMask, :);
    
    assert(length(unique(visTimes)) == length(visTimes), ...
        'Repeated visits while padding data!');
    assert(all(samplingTimes > 0), ...
        'Data is not sorted or there are duplicate visits!');
    
    if ~all(abs(samplingTimes) == res)
        idealVisits     = visTimes(1):res:visTimes(end);
        nIdealVisits    = length(idealVisits);
        isObserved      = ismember(idealVisits, visTimes);
        obsIdx          = find(isObserved);
        
        newData         = num2cell(nan(nIdealVisits, size(patData, 2)));
        newPhi          = nan(nIdealVisits, size(patPhi, 2));
        
        newData(isObserved, :)      = patData;
        newPhi(isObserved, :)       = patPhi;
        
        % linearly interpolate some columns
        for ii = 1:length(linDataColIdx)
            colIdx          = linDataColIdx(ii);
            col             = cat(1, newData{:, colIdx});
            newCol          = interp1(idealVisits(isObserved), ...
                col(isObserved), idealVisits, 'linear');
            newData(~isObserved, colIdx)    = num2cell(newCol(~isObserved));
        end
        
        % for others like RID, apply zero order hold
        for ii = 1:length(holdDataColIdx)
            colIdx          = holdDataColIdx(ii);
            col             = newData(:, colIdx);
            
            if isnumeric(col{obsIdx(1)})
                col         = cat(1, col{:});
                newCol      = interp1(idealVisits(isObserved), ...
                    col(isObserved), idealVisits, 'previous');
                newData(~isObserved, colIdx)    = ...
                    num2cell(newCol(~isObserved));
            else
                fillVal     = repmat(col(obsIdx(1)), sum(~isObserved), 1);
                newData(~isObserved, colIdx)    = fillVal;
            end
        end

        % interpolate diagnosis only if patient has the same dx before
        % and after the period of missingness
        for uu=1:(length(obsIdx)-1)
            beginIdx        = obsIdx(uu);
            endIdx          = obsIdx(uu + 1);
            if (endIdx - beginIdx) > 1 && ...
                    (newData{beginIdx, diagnosisIdx} == ...
                    newData{endIdx, diagnosisIdx})
                newData(beginIdx+1:endIdx-1, diagnosisIdx)  = ...
                    repmat(newData(beginIdx, diagnosisIdx), ...
                    endIdx-beginIdx-1, 1);
            end
        end
        patData         = newData;
        patPhi          = newPhi;
    end
    
    paddedData{p}       = patData;
    paddedPhi{p}        = patPhi;
end

newData     = cell2table(cat(1, paddedData{:}), ...
    'VariableNames', dataColNames);
newPhi      = cat(1, paddedPhi{:});

% newPhi(:, strcmpi(allData.phiNames, 'AGE'))     = newData.AGE;
% newPhi(:, strcmpi(allData.phiNames, 'ApoE-Profile'))    = ...
%     getGeneticSummary(newData);
% newPhi(:, strcmpi(allData.phiNames, 'EDUCATION'))       = newData.PTEDUCAT;

allData.data    = newData;
allData.phi     = newPhi;

end