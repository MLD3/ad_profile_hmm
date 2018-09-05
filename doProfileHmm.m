% load profile-hmm-toy-data.mat;
% trajectories    = observations;

load profile-hmm-data.mat;

%% pick out the clinical score used to intialize the data; we use ADAS-Cog
adasIdx             = strcmpi('ADAS', allData.phiNames);
score               = cellfun(@(seq)seq(:, adasIdx), trajectories, ...
    'UniformOutput', false);

%% Discretize the data
rawTrajectories     = trajectories;

[rawTrajectories, dx, yl, yu, followup, score]  = ...
    padHmmTrajectories(rawTrajectories, dx, yl, yu, followup, score, 6);

stackedPhi          = cat(1, rawTrajectories{:});
counts              = cellfun(@length, yl);

nDiscreteValues     = 10;
[discretePhi, quantileLims]     = getProjectedFeatures(stackedPhi, ...
    'quantile', 'nQuantiles', nDiscreteValues, ...
    'binaryOutput', false, 'isMissing', isnan(stackedPhi), ...
    'removeRedundant', false, 'missingInd', false);
    
trajectories        = mat2cell(discretePhi, counts);
data                = cellfun(@transpose, trajectories, ...
    'UniformOutput', false);

%% Setup the experiment parameters
kRange          = [4:2:30];
nExpRange       = [0 1 2 3 4 5 6];
initMethod      = {'score'};

[kGrid, nExpGrid, methodGrid]   = ndgrid(kRange, nExpRange, ...
    1:length(initMethod));
methodGrid              = initMethod(methodGrid);

nModels                 = numel(kGrid);

%% Save both the baseline model, as well as the trained HMM model
models                  = cell(size(kGrid));
paths                   = cell(size(kGrid));

hoLikelihood            = zeros(length(trajectories), numel(kGrid));
likelihood              = zeros(length(trajectories), numel(kGrid));

parfor exp=1:nModels
    nStatesOpt      = kGrid(exp);
    nExpOpt         = nExpGrid(exp);
    methodOpt       = methodGrid{exp};

    loglik          = zeros(size(trajectories));
    
    for f=1:hofolds.NumTestSets
        trainIdx        = hofolds.training(:, f);
        testIdx         = hofolds.test(:, f);
        
        initModel       = getInitModel(trajectories(trainIdx), ...
            nStatesOpt, 'discrete', 'method', methodOpt, ...
            'nDiscreteValues', nDiscreteValues, 'score', ...
            score(trainIdx), 'nExpansions', nExpOpt);

        [model, loglikHist]     = hmmFitEm(data(trainIdx), ...
            initModel.nstates, 'discrete', 'model', initModel, ...
            'verbose', true);
        
        loglik(testIdx)     = hmmLogprob(model, data(testIdx));
    end
    
    hoLikelihood(:, exp)        = loglik;

    initModel   = getInitModel(trajectories, nStatesOpt, 'discrete', ...
        'method', methodOpt, 'score', score, 'nExpansions', nExpOpt, ...
        'nDiscreteValues', nDiscreteValues);
    
    model       = hmmFitEm(data, initModel.nstates, 'discrete', ...
        'model', initModel, 'verbose', true);
    likelihood(:, exp)      = hmmLogprob(model, data);
    
    path            = cellfun(@(seq)hmmMap(model, seq), data, ...
        'UniformOutput', false);

    models{exp}             = model;
    paths{exp}              = path;
end

save('profile-hmm-es-discrete');