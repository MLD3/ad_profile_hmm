% load profile-hmm-toy-data.mat;
% trajectories    = observations;

load profile-hmm-data.mat;

%% pick out the clinical score used to intialize the data; we use ADAS-Cog
adasIdx             = strcmpi('ADAS', allData.phiNames);
score               = cellfun(@(seq)seq(:, adasIdx), trajectories, ...
    'UniformOutput', false);

%% Discretize the data
rawTrajectories     = trajectories;
[rawTrajectoriesPad, dxPad, ylPad, yuPad, followupPad, scorePad]  = ...
    padHmmTrajectories(trajectories, dx, yl, yu, followup, score, 6);

rawPhi              = cat(1, trajectories{:});
rawPhiPad           = cat(1, rawTrajectoriesPad{:});
counts              = cellfun(@(seq)size(seq, 1), rawTrajectories);
countsPad           = cellfun(@(seq)size(seq, 1), rawTrajectoriesPad);

nDiscreteValues     = 10;
skipDiscrete        = false(length(allData.phiNames), 1);
skipDiscrete(strcmpi(allData.phiNames, 'ApoE-Profile'))     = true;

%% Ignore missing data, pretend it is all uniformly samples
discretePhi         = getProjectedFeatures(rawPhi, ...
    'quantile', 'nQuantiles', nDiscreteValues, ...
    'binaryOutput', false, 'isMissing', isnan(rawPhi), ...
    'removeRedundant', false, 'skipDiscrete', skipDiscrete);

trajectories        = mat2cell(discretePhi, counts);
data                = cellfun(@transpose, trajectories, ...
    'UniformOutput', false);

%% Explicitly model missing data with a discrete indicator (MCAR)
discretePhi         = getProjectedFeatures(rawPhiPad, ...
    'quantile', 'nQuantiles', nDiscreteValues, ...
    'binaryOutput', false, 'isMissing', isnan(rawPhiPad), ...
    'removeRedundant', false, 'skipDiscrete', skipDiscrete);

trajectoriesPad     = mat2cell(discretePhi, countsPad);
dataPad             = cellfun(@transpose, trajectoriesPad, ...
    'UniformOutput', false);

%% Assume missingness at random (MAR) and marginalize missed observations
discretePhi         = getProjectedFeatures(rawPhiPad, ...
    'quantile', 'nQuantiles', nDiscreteValues, ...
    'binaryOutput', false, 'isMissing', isnan(rawPhiPad), ...
    'removeRedundant', false, 'missingInd', false, ...
    'skipDiscrete', skipDiscrete);

trajectoriesNan     = mat2cell(discretePhi, countsPad);
dataNan             = cellfun(@transpose, trajectoriesNan, ...
    'UniformOutput', false);

counts      = countsPad;

%% Setup the experiment parameters
kGrid           = 12:2:30;
nModels         = length(kGrid);

nBaselineModels         = numel(kGrid);

missingLoglik           = zeros(length(data), nBaselineModels);
paddedLoglik            = zeros(length(data), nBaselineModels);
nanLoglik               = zeros(length(data), nBaselineModels);

missingLoglikZ          = zeros(length(data), nBaselineModels);
paddedLoglikZ           = zeros(length(data), nBaselineModels);
nanLoglikZ              = zeros(length(data), nBaselineModels);

missingModels           = cell(nBaselineModels, 1);
paddedModels            = cell(nBaselineModels, 1);
nanModels               = cell(nBaselineModels, 1);

missingPaths            = cell(nBaselineModels, 1);
paddedPaths             = cell(nBaselineModels, 1);
nanPaths                = cell(nBaselineModels, 1);

for exp=1:nBaselineModels
    nStatesOpt      = kGrid(exp);
    
    missLL          = zeros(length(trajectories), 1);
    padLL           = zeros(length(trajectories), 1);
    nanLL           = zeros(length(trajectories), 1);
    
    missLLZ         = zeros(length(trajectories), 1);
    padLLZ          = zeros(length(trajectories), 1);
    nanLLZ          = zeros(length(trajectories), 1);
    
    for f=1:hofolds.NumTestSets
        trainIdx        = hofolds.training(:, f);
        testIdx         = hofolds.test(:, f);
        
        initModel       = getInitModel(trajectories(trainIdx), ...
            nStatesOpt, 'discrete', 'method', 'score', ...
            'nDiscreteValues', nDiscreteValues, 'score', ...
            score(trainIdx), 'nExpansions', 0, 'dx', dx(trainIdx), ...
            'sortModel', true);
        
        model           = hmmFitEm(data(trainIdx), ...
            initModel.nstates, 'discrete', 'model', initModel, ...
            'verbose', true);
        
        trainLL             = hmmLogprob(model, dataNan(trainIdx));
        missLL(testIdx)     = hmmLogprob(model, dataNan(testIdx));
        missLLZ(testIdx)    = getLoglikZscore(missLL, counts, trainIdx, ...
            'trainLL', trainLL, 'testIdx', testIdx);
        
        initModel       = getInitModel(trajectoriesPad(trainIdx), ...
            nStatesOpt, 'discrete', 'method', 'score', ...
            'nDiscreteValues', nDiscreteValues, 'score', ...
            scorePad(trainIdx), 'nExpansions', 0, ...
            'modelMissingness', true, 'dx', dxPad(trainIdx), ...
            'sortModel', true);
        
        [model, ~]      = hmmFitEm(dataPad(trainIdx), ...
            initModel.nstates, 'discrete', 'model', initModel, ...
            'verbose', true);
        
        trainLL             = hmmLogprob(model, dataNan(trainIdx));
        padLL(testIdx)      = hmmLogprob(model, dataNan(testIdx));
        padLLZ(testIdx)     = getLoglikZscore(padLL, counts, trainIdx, ...
            'trainLL', trainLL, 'testIdx', testIdx);
        
        initModel       = getInitModel(trajectoriesNan(trainIdx), ...
            nStatesOpt, 'discrete', 'method', 'score', ...
            'nDiscreteValues', nDiscreteValues, 'score', ...
            scorePad(trainIdx), 'nExpansions', 0, 'dx', dxPad(trainIdx), ...
            'sortModel', true);
        
        [model, ~]      = hmmFitEm(dataNan(trainIdx), ...
            initModel.nstates, 'discrete', 'model', initModel, ...
            'verbose', true);
        
        trainLL             = hmmLogprob(model, dataNan(trainIdx));
        nanLL(testIdx)      = hmmLogprob(model, dataNan(testIdx));
        nanLLZ(testIdx)     = getLoglikZscore(nanLL, counts, trainIdx, ...
            'trainLL', trainLL, 'testIdx', testIdx);
    end
    
    missingLoglik(:, exp)       = missLL;
    paddedLoglik(:, exp)        = padLL;
    nanLoglik(:, exp)           = nanLL;
    
    missingLoglikZ(:, exp)      = missLLZ;
    paddedLoglikZ(:, exp)       = padLLZ;
    nanLoglikZ(:, exp)          = nanLLZ;

    initModel       = getInitModel(trajectories, nStatesOpt, ...
        'discrete', 'method', 'score', 'nExpansions', 0, ...
        'nDiscreteValues', nDiscreteValues, 'score', score, ...
        'dx', dx, 'sortModel', true);
    missingModel    = hmmFitEm(data, initModel.nstates, 'discrete', ...
        'model', initModel, 'verbose', true);
    
    initModel       = getInitModel(trajectoriesPad, nStatesOpt, ...
        'discrete', 'method', 'score', 'nExpansions', 0, ...
        'nDiscreteValues', nDiscreteValues, 'score', scorePad, ...
        'dx', dxPad, 'sortModel', true);
    paddedModel     = hmmFitEm(dataPad, initModel.nstates, 'discrete', ...
        'model', initModel, 'verbose', true);
    
    initModel       = getInitModel(trajectoriesNan, nStatesOpt, ...
        'discrete', 'method', 'score', 'nExpansions', 0, ...
        'nDiscreteValues', nDiscreteValues, 'score', scorePad, ...
        'dx', dxPad, 'sortModel', true);
    nanModel        = hmmFitEm(dataNan, initModel.nstates, 'discrete', ...
        'model', initModel, 'verbose', true);
    
    missingPaths{exp}   = cellfun(@(seq)hmmMap(missingModel, seq), ...
        dataNan, 'UniformOutput', false);
    paddedPaths{exp}    = cellfun(@(seq)hmmMap(paddedModel, seq), ...
        dataNan, 'UniformOutput', false);
    nanPaths{exp}       = cellfun(@(seq)hmmMap(nanModel, seq), ...
        dataNan, 'UniformOutput', false);
    
    missingModels{exp}      = missingModel;
    paddedModels{exp}       = paddedModel;
    nanModels{exp}          = nanModel;
end

save('profile-hmm-missingness');