function [allData, hasL, hasU, censored] = prepMultiModalData(varargin)

[   groups                                              , ...
    scale                                               , ...
    summaryData                                         , ...
    useUncAd                                            , ...
    useAnyAd                                            , ...
    useNl                                               , ...
    removeZeroFollowup                                  , ...
    baseGroup                                           , ...
    minVisits                                           , ...
    leftCensored                                        , ...
    rightCensored           ] = process_options(varargin, ...
    'groups'                , {{'MRI', 'DEMOGRAPHIC', 'GENETIC'}, ...
    {'FDG'}, {'CSF'}, {'AMYLOID'}}                      , ...
    'scale'                 , true                      , ...
    'summaryData'           , false                     , ...
    'useUncAd'              , false                     , ...
    'useAnyAd'              , true                      , ...
    'useNl'                 , false                     , ...
    'removeZeroFollowup'    , true                      , ...
    'baseGroup'             , []                        , ...
    'minVisits'             , 1                         , ...
    'leftCensored'          , true                      , ...
    'rightCensored'         , true                      );

NL              = 1;
MCI             = 2;
AD              = 3;

if summaryData
    modalities      = num2cell(cat(2, groups{:}));
    
    [phiGrouped, data, profileIdx, phiNames, phiType, uniqueGroups] = ...
        getMultiModalData('minVisits', 1, 'groups', modalities, ...
        'scale', false);
    
    for g=1:length(modalities)
        if strcmpi(modalities{g}{1}, 'MRI')
            [phiGrouped{g}, ~, phiNames{g}]    = getMriSummary(data, ...
                'scale', false, 'imputeNan', false, 'age', false, ...
                'normIcv', false);
        elseif strcmpi(modalities{g}{1}, 'FDG')
            [phiGrouped{g}, ~, phiNames{g}]    = getPetSummary(data, ...
                'scale', false, 'imputeNan', false, 'age', false);
        elseif strcmpi(modalities{g}{1}, 'CLINICAL')
            [phiGrouped{g}, ~, phiNames{g}]    = getClinicalSummary(data, ...
                'scale', false, 'imputeNan', false, 'age', false);
        elseif strcmpi(modalities{g}{1}, 'AMYLOID')
            phiGrouped{g}       = data{:, 'SUMMARYSUVR_WHOLECEREBNORM'};
            phiNames{g}         = {'AMYLOID'};
        elseif strcmpi(modalities{g}{1}, 'GENETIC')
            [phiGrouped{g}, ~, phiNames{g}]     = getGeneticSummary(...
                data);
        elseif strcmpi(modalities{g}{1}, 'DEMOGRAPHIC')
            phiGrouped{g}       = data{:, {'AGE', 'PTEDUCAT'}};
            phiNames{g}         = {'AGE', 'EDUCATION'};
        end
    end
    
    isinGroup       = zeros(size(data, 1), length(groups));
    phiGroupedNew   = cell(length(groups), 1);
    phiNamesNew     = cell(length(groups), 1);
    
    for g=1:length(groups)
        groupModalities     = groups{g};
        modalityMask        = false(length(modalities), 1);
        
        for m=1:length(groupModalities)
            modalityMask    = modalityMask | ...
                strcmpi(cat(1, modalities{:}), groupModalities{m});
        end
        
        matches             = find(all(uniqueGroups(:, modalityMask), 2));
        isinGroup(:, g)     = ismember(profileIdx, matches);
        
        phiGroupedNew{g}    = cat(2, phiGrouped{modalityMask});
        phiNamesNew{g}      = cat(2, phiNames{modalityMask});
    end
    
    [uniqueGroups, ~, profileIdx]   = unique(isinGroup, 'rows');
    phiGrouped      = phiGroupedNew;
    phiNames        = phiNamesNew;
else
    [phiGrouped, data, profileIdx, phiNames, uniqueGroups] = ...
        getMultiModalData('minVisits', 1, 'groups', groups, 'scale', false);
end

d               = cellfun(@(seq)size(seq, 2), phiGrouped);

groupIdx        = arrayfun(@(g)ones(1, d(g))*g, 1:length(d), ...
    'UniformOutput', false);
groupIdx        = cat(2, groupIdx{:});

phi             = cat(2, phiGrouped{:});
phiNames        = cat(2, phiNames{:});

allData.phiNames        = phiNames;
allData.phiType         = phiType;
allData.phi             = phi;
allData.group           = groupIdx;

index           = true(size(data.DX));

if ~useNl
    index       = index & (~(data.DX == NL));
end

trivial                 = (data.CONVTIME_L == -999 & ...
    data.CONVTIME_U == 999);

index                   = index & ~trivial;

hasL                    = data.CONVTIME_L > -999;
hasU                    = data.CONVTIME_U < 999;

if ~leftCensored
    index               = index & hasL;
end

if ~rightCensored
    index               = index & hasU;
end

if ~useUncAd
    index               = index & ~(data.DX == AD & hasL);
end

if ~useAnyAd
    index               = index & (data.DX == MCI);
end

if removeZeroFollowup
    index               = index & (data.FOLLOWUP > 0);
end

if ~isempty(baseGroup)
    index               = index & ~any(isnan(phi(:, ...
        ismember(groupIdx, baseGroup))), 2);
end

[patients, ~, patientId]    = unique(data.RID);
counts                      = histc(data.RID, patients);
counts                      = counts(patientId);
index                       = index & (counts >= minVisits);

if scale
    phi                 = zScale(phi);
end

allData.phi             = phi(index, :);
allData.data            = data(index, :);
allData.profile         = profileIdx(index);
allData.profileMakeup   = uniqueGroups;

data                    = allData.data;
convtimeL               = data.CONVTIME_L;
convtimeU               = data.CONVTIME_U;

hasL                    = convtimeL > -999;
hasU                    = convtimeU < 999;

% treat uncensored AD as censored because we don't want to make predictions
% for those
censored                = ~hasL | ~hasU;
    
end