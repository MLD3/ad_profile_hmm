% load profile-hmm-toy-data.mat;
% trajectories    = observations;

load profile-hmm-data.mat;

%% pick out the clinical score used to intialize the data; we use ADAS-Cog
adasIdx             = strcmpi('ADAS', allData.phiNames);
score               = cellfun(@(seq)seq(:, adasIdx), trajectories, ...
    'UniformOutput', false);

%% Discretize the data
rawTrajectories     = trajectories;

rawTrajectoriesUnc  = rawTrajectories(~censored);
scoreUnc            = score(~censored);
ylUnc               = yl(~censored);
yuUnc               = yu(~censored);
dxUnc               = dx(~censored);
followupUnc         = followup(~censored);

[rawTrajectoriesUnc, dxUnc, ylUnc, yuUnc, followupUnc, scoreUnc]    = ...
    padHmmTrajectories(rawTrajectoriesUnc, dxUnc, ylUnc, yuUnc, ...
    followupUnc, scoreUnc, 6);

stackedPhi          = cat(1, rawTrajectoriesUnc{:});
counts              = cellfun(@(seq)size(seq, 1), rawTrajectoriesUnc);

nDiscreteValues     = 10;
[discretePhi, quantileLims]     = getProjectedFeatures(stackedPhi, ...
    'quantile', 'nQuantiles', nDiscreteValues, ...
    'binaryOutput', false, 'isMissing', isnan(stackedPhi), ...
    'removeRedundant', false, 'missingInd', false);

trajectoriesUnc     = mat2cell(discretePhi, counts);
dataUnc             = cellfun(@transpose, trajectoriesUnc, ...
    'UniformOutput', false);

%% Setup the experiment parameters
kGrid           = [3 4 5 6 7 8 9 10 12 14 16];
nModels         = length(kGrid);

nBaselineModels         = numel(kGrid);

convtimeLoglik          = zeros(length(dataUnc), nBaselineModels);
tempConvLoglik          = zeros(length(dataUnc), nBaselineModels);
hmmLoglik               = zeros(length(dataUnc), nBaselineModels);

convtimeModels          = cell(nBaselineModels, 1);
tempConvModels          = cell(nBaselineModels, 1);
hmmModels               = cell(nBaselineModels, 1);

parfor exp=1:nBaselineModels
    nStatesOpt      = kGrid(exp);
    
    convLL          = zeros(length(trajectoriesUnc), 1);
    tempConvLL      = zeros(length(trajectoriesUnc), 1);
    hmmLL           = zeros(length(trajectoriesUnc), 1);
    
    for f=1:hofolds.NumTestSets
        trainIdx        = hofoldsUnc.training(:, f);
        testIdx         = hofoldsUnc.test(:, f);
        
        convtimeModel   = getConvtimeModel(trajectoriesUnc(trainIdx), ...
            nStatesOpt, 'discrete', ylUnc(trainIdx), yuUnc(trainIdx), ...
            'nDiscreteValues', nDiscreteValues);
        convLL(testIdx)     = convtimeModelLogProb(convtimeModel, ...
            trajectoriesUnc(testIdx), ylUnc(testIdx), yuUnc(testIdx));
        
        convModel       = getInitModel(trajectoriesUnc(trainIdx), ...
            nStatesOpt, 'discrete', 'method', 'convtime', ...
            'nDiscreteValues', nDiscreteValues, 'yl', ylUnc(trainIdx), ...
            'yu', yuUnc(trainIdx), 'nExpansions', 0, ...
            'modelMissingness', false);
        tempConvLL(testIdx)     = hmmProfileLogprob(convModel, ...
            dataUnc(testIdx));
        
        initModel       = getInitModel(trajectoriesUnc(trainIdx), ...
            nStatesOpt, 'discrete', 'method', 'score', ...
            'nDiscreteValues', nDiscreteValues, 'score', ...
            scoreUnc(trainIdx), 'nExpansions', 0, 'modelMissingness', ...
            false);
        
        [model, loglikHist]     = hmmProfileFitEm(dataUnc(trainIdx), ...
            'condDiscreteProd', 'model', initModel, ...
            'EMargs', {'verbose', true});
        
        hmmLL(testIdx)      = hmmProfileLogprob(model, dataUnc(testIdx));
    end
    
    convtimeLoglik(:, exp)      = convLL;
    hmmLoglik(:, exp)           = hmmLL;
    tempConvLoglik(:, exp)      = tempConvLL;
    
    convtimeModel   = getConvtimeModel(trajectoriesUnc, nStatesOpt, ...
        'discrete', ylUnc, yuUnc, 'nDiscreteValues', nDiscreteValues);
    
    convModel       = getInitModel(trajectoriesUnc, ...
        nStatesOpt, 'discrete', 'method', 'convtime', ...
        'nDiscreteValues', nDiscreteValues, 'yl', ylUnc, ...
        'yu', yuUnc, 'nExpansions', 0, 'modelMissingness', false);
    
    initModel       = getInitModel(trajectoriesUnc, nStatesOpt, ...
        'discrete', 'method', 'score', 'nExpansions', 0, ...
        'nDiscreteValues', nDiscreteValues, 'score', scoreUnc, ...
        'modelMissingness', false);
    hmmModel        = hmmProfileFitEm(dataUnc, 'condDiscreteProd', ...
        'model', initModel, 'EMargs', {'verbose', true});
    
    convtimeModels{exp}     = convtimeModel;
    tempConvModels{exp}     = convModel;
    hmmModels{exp}          = hmmModel;
end

save('profile-hmm-vs-baseline');