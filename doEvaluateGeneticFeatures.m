load profile-hmm-data.mat;

%% pick out the clinical score used to intialize the data; we use ADAS-Cog
scoreIdx            = strcmpi('ADAS', allData.phiNames);
score               = cellfun(@(seq)seq(:, scoreIdx), trajectories, ...
    'UniformOutput', false);

meanScore           = nanmean(cat(1, score{:}));

for ii=1:length(score)
    patScore        = score{ii};
    patFollow       = followup{ii};
    
    obsScore        = ~isnan(patScore);
    
    if sum(obsScore) >= 2
        score{ii}   = interp1(patFollow(obsScore), patScore(obsScore), ...
            patFollow, 'linear', 'extrap');
    elseif sum(obsScore) == 1
        score{ii}   = ones(size(obsScore))*patScore(obsScore);
    else
        score{ii}   = ones(size(obsScore))*meanScore;
    end
end

%% Discretize the data
rawTrajectories     = trajectories;

rawPhi              = cat(1, rawTrajectories{:});
counts              = allData.counts;

nDiscreteValues     = 10;
skipDiscrete        = false(length(allData.phiNames), 1);
skipDiscrete(strcmpi(allData.phiNames, 'ApoE-Profile'))     = true;
% skipDiscrete(strcmpi(allData.phiNames, 'EDUCATION'))     = true;

%% Assume missingness at random (MAR) and marginalize missed observations
discretePhi         = getProjectedFeatures(rawPhi, ...
    'quantile', 'nQuantiles', nDiscreteValues, ...
    'binaryOutput', false, 'isMissing', isnan(rawPhi), ...
    'removeRedundant', false, 'missingInd', false, ...
    'skipDiscrete', skipDiscrete);

geneIdx             = strcmpi(allData.phiNames, 'ApoE-Profile');
demogIdx            = strcmpi(allData.phiNames, 'AGE') | ...
    strcmpi(allData.phiNames, 'EDUCATION');

trajectories        = mat2cell(discretePhi(:, ~geneIdx & ~demogIdx), counts);
trajectoriesDemog   = mat2cell(discretePhi(:, ~geneIdx), counts);
trajectoriesGene    = mat2cell(discretePhi(:, ~demogIdx), counts);
trajectoriesBoth    = mat2cell(discretePhi, counts);

data                = cellfun(@transpose, trajectories, ...
    'UniformOutput', false);
dataDemog           = cellfun(@transpose, trajectoriesDemog, ...
    'UniformOutput', false);
dataGene            = cellfun(@transpose, trajectoriesGene, ...
    'UniformOutput', false);
dataBoth            = cellfun(@transpose, trajectoriesBoth, ...
    'UniformOutput', false);

%% Setup the experiment parameters
kGrid               = 4:2:30;

trajOptions         = {trajectories, trajectoriesDemog, trajectoriesGene, ...
    trajectoriesBoth};
dataOptions         = {data, dataDemog, dataGene, dataBoth};

[kOpt, trajOpt]     = ndgrid(kGrid, 1:length(trajOptions));
dataOpt             = dataOptions(trajOpt);
trajOpt             = trajOptions(trajOpt);

nModels             = numel(kOpt);

loglik              = zeros(length(data), nModels);
loglikZ             = zeros(length(data), nModels);

models              = cell(nModels, 1);

parfor exp=1:nModels
    nStatesOpt      = kOpt(exp);
    traj            = trajOpt{exp};
    dat             = dataOpt{exp};
    
    ll              = zeros(length(traj), 1);
    llZ             = zeros(length(traj), 1);
    
    for f=1:hofolds.NumTestSets
        trainIdx        = hofolds.training(:, f);
        testIdx         = hofolds.test(:, f);
            
        initModel       = getInitModel(traj(trainIdx), ...
            nStatesOpt, 'discrete', 'method', 'score', ...
            'nDiscreteValues', nDiscreteValues, 'score', ...
            score(trainIdx), 'nExpansions', 0, 'dx', dx(trainIdx), ...
            'sortModel', false);
            
        [model, ~]      = hmmFitEm(dat(trainIdx), ...
            initModel.nstates, 'discrete', 'model', initModel, ...
            'verbose', true);
            
        trainLL             = hmmLogprob(model, dat(trainIdx));
        ll(testIdx)         = hmmLogprob(model, dat(testIdx));
        llZ(testIdx)        = getLoglikZscore(ll, counts, ...
            trainIdx, 'trainLL', trainLL, 'testIdx', testIdx);
    end
    
    loglik(:, exp)      = ll;
    loglikZ(:, exp)     = llZ;

    initModel       = getInitModel(traj, nStatesOpt, ...
        'discrete', 'method', 'score', 'nExpansions', 0, ...
        'nDiscreteValues', nDiscreteValues, 'score', score, ...
        'dx', dx, 'sortModel', false);
    [model, llHistS]    = hmmFitEm(dat, ...
        initModel.nstates, 'discrete', 'model', initModel, ...
        'verbose', true);
        
    models{exp}     = model;
end

loglik          = reshape(loglik, length(trajectories), length(kGrid), []);
loglikZ         = reshape(loglikZ, length(trajectories), length(kGrid), []);
models          = reshape(models, length(kGrid), []);

save('profile-hmm-feature-comparison');