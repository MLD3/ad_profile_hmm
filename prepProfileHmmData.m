% what kinds of data do we want?
sources     = {{'MRI'}, {'FDG'}, {'CSF'}, {'AMYLOID'}, {'CLINICAL'}, ...
    {'GENETIC'}, {'DEMOGRAPHIC'}};

[allData, hasL, hasU, censored]     = prepMultiModalData(...
    'groups', sources, 'useAnyAd', true, 'useUncAd', true, ...
    'leftCensored', true, 'rightCensored', true, 'scale', false, ...
    'summaryData', true, 'removeZeroFollowup', false, 'useNl', true);

allData     = prepTemporalData(allData, 'minCount', 2, ...
    'noMissVar', {}, 'resample', true, 'samplingRate', 6);

% fillup the binary missingness matrix that tells us about the modalities
% available for each visit
isMissing               = isnan(allData.phi);
allData.phiUnscaled     = allData.phi;

notClinical             = allData.group ~= 5;
allData.phi             = allData.phiUnscaled;
allData.phi(:, notClinical)     = zScale(allData.phi(:, notClinical));

allData.phi(isMissing)  = nan;

data            = allData.data;
trajIdx         = allData.trajIdx;
nTraj           = length(trajIdx);

trajectories    = colvec(arrayfun(@(ii)allData.phiUnscaled(trajIdx{ii}, :), ...
    1:nTraj, 'UniformOutput', false));
yl              = colvec(arrayfun(@(ii)data.CONVTIME_L(trajIdx{ii}, :), ...
    1:nTraj, 'UniformOutput', false));
yu              = colvec(arrayfun(@(ii)data.CONVTIME_U(trajIdx{ii}, :), ...
    1:nTraj, 'UniformOutput', false));
followup        = colvec(arrayfun(@(ii)data.FOLLOWUP(trajIdx{ii}, :), ...
    1:nTraj, 'UniformOutput', false));
dx              = colvec(arrayfun(@(ii)data.DX(trajIdx{ii}, :), ...
    1:nTraj, 'UniformOutput', false));
age             = colvec(arrayfun(@(ii)data.AGE(trajIdx{ii}, :), ...
    1:nTraj, 'UniformOutput', false));

censored        = cellfun(@(lower, upper)any(lower == -999) || ...
    any(upper == 999), yl, yu);

counts          = cellfun(@length, yl);
firstYl         = cellfun(@(seq)seq(1), yl);
firstYu         = cellfun(@(seq)seq(1), yu);

convtime                    = firstYu;
convtime(convtime == 999)   = firstYl(convtime == 999);
rid                         = data.RID;
[patients, ~, ridGroup]     = unique(rid);

isVisitMissing  = all(isnan(allData.phiUnscaled), 2);

convtimeQ       = discretize(convtime, ...
    quantile(convtime, linspace(0, 1, 6)));
countQ          = discretize(counts, quantile(counts, linspace(0, 1, 6)));

[~, ~, groups]  = unique([convtimeQ countQ], 'rows');

%% Setup hofolds that will be used to calculate held-out LL for the HMM
hofolds         = get_folds('lkpo'      , ...
    'groups'    , groups                , ...
    'patients'  , patients              , ...
    'rid'       , patients              , ...
    'nRep'      , 1                     , ...
    'nFolds'    , 10                    , ...
    'getSparse' , false                 );

%% Leave-One-patient-out hofolds for the time-series
% hofSeries       = struct('NumTestSets', length(patients), ...
%     'test', logical(eye(length(patients))), 'patientSet', ...
%     ones(1, length(patients)), 'type', 'lopo');
% hofSeries.training  = ~(hofSeries.test);

hofSeries       = get_folds('lkpo'      , ...
    'groups'    , groups                , ...
    'patients'  , patients              , ...
    'rid'       , patients              , ...
    'nRep'      , 1                     , ...
    'nFolds'    , 10                    , ...
    'getSparse' , false                 );

nHof            = hofSeries.NumTestSets;

%% Leave one (or k) patient(s) out, but now set them up on an instance level
hofInst         = struct(hofSeries);
training        = false(length(rid), nHof);
test            = false(length(rid), nHof);

for f=1:nHof
    training(ismember(rid, patients(hofSeries.training(:, f))), f)  = true;
    test(ismember(rid, patients(hofSeries.test(:, f))), f)          = true;
end

hofInst.training    = training;
hofInst.test        = test;

assert(hofInst.NumTestSets == hofSeries.NumTestSets, ...
    'Mismatch in folds between time-series and data instances!');

%% setup individual hold-out folds for each conversion time horizon
durationGrid        = [12 24 36 48 60 72];

nInst               = size(hofInst.training, 1);
nDurations          = length(durationGrid);

hofClf              = cell(length(durationGrid), 1);
cvfClf              = cell(length(durationGrid), 1);
cvfHandleClf        = cell(length(durationGrid), 1);

labels              = zeros(nInst, nDurations);
censored            = zeros(nInst, nDurations);

for d=1:nDurations
    temp            = getConversionWithinDuration(data, durationGrid(d));
    labels(:, d)    = temp.CONV;
    censored(:, d)  = temp.CENSORED;
    
    classificationMask      = ~censored(:, d);
    
    patientMask     = ismember(patients, rid(classificationMask));
    durPatients     = patients(patientMask);
    durRid          = rid(classificationMask);
    
    durCounts       = histc(durRid, durPatients);
    patLabels       = mat2cell(labels(classificationMask, d), durCounts);
    groups          = cellfun(@(seq)sum(seq == 1), patLabels);
    
    hof             = get_folds('lkpo'      , ...
        'patients'  , durPatients           , ...
        'groups'    , groups                , ...
        'rid'       , durRid                , ...
        'nFolds'    , 10                    , ...
        'nRep'      , 1                     , ...
        'getSparse' , true                  );
    nHof            = hof.NumTestSets;
    
    cvfHandle       = @(foldNum, hof)get_cvfolds(foldNum, hof, ...
        durRid, durPatients, groups, 'lkpo', 5, 1, 0.2);
    
    hofClf{d}           = hof;
    cvfClf{d}           = repmat({{}}, hof.NumTestSets, 1);
    cvfHandleClf{d}     = cvfHandle;
end

%% Setup folds for the regression problem

phiToPredict        = strcmpi(allData.phiNames, 'ADAS');
target              = allData.phi(:, phiToPredict);
regressionMask      = ~isnan(target);
target              = target(regressionMask);

patientMask         = ismember(patients, rid(regressionMask));
regPatients         = patients(patientMask);
regRid              = rid(regressionMask);

regCounts           = histc(regRid, regPatients);
patTarget           = mat2cell(target, regCounts);
groups              = cellfun(@mean, patTarget);
groups              = discretize(groups, quantile(groups, ...
    linspace(0, 1, 6)));
    
hofReg              = get_folds('lkpo'      , ...
    'patients'      , regPatients           , ...
    'groups'        , groups                , ...
    'rid'           , regRid                , ...
    'nFolds'        , 10                    , ...
    'nRep'          , 1                     , ...
    'getSparse'     , true                  );

cvfReg              = repmat({{}}, hofReg.NumTestSets, 1);
    
cvfHandleReg        = @(foldNum, hof)get_cvfolds(foldNum, hof, ...
        regRid, regPatients, groups, 'lkpo', 5, 1, 0.2);

save('profile-hmm-data.mat');