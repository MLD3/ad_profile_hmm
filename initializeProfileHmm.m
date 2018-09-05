function initModel = initializeProfileHmm(data, ...
    varargin)
%% Initialize a HMM Model for global sequence alignment
% data          - n x 1 cell, each cell i is a ni x d matrix

[   nStates                                             , ...
    cpdType                                             , ...
    nDiscreteValues                                     , ...
    yl                                                  , ...
    yu                                                  , ...
    initMethod                                          , ...
    score                   ] = process_options(varargin, ...
    'nStates'               , 20                        , ...
    'cpdType'               , 'discrete'                , ...
    'nDiscreteValues'       , 10                        , ...
    'yl'                    , {}                        , ...
    'yu'                    , {}                        , ...
    'initMethod'            , 'mixture'                 , ...
    'score'                 , {}                        );

stackedPhi      = cat(1, data{:});
counts          = cellfun(@(seq)size(seq, 1), data);

if ~isempty(score)
    stackedScore    = cat(1, score{:});
    
    binLimits       = quantile(stackedScore, linspace(0, 1, nStates+1));
    stackedStates   = discretize(stackedScore, binLimits);
    
    states          = mat2cell(stackedStates, counts);
else
    if strcmpi(initMethod, 'mixture')
        stackedStates   = getInitStates(stackedPhi, nStates, cpdType);
    elseif strcmpi(initMethod, 'convtime')
        stackedYl       = cat(1, yl{:});
        stackedYu       = cat(1, yu{:});
        
        lCensored       = stackedYl == -999;
        rCensored       = stackedYu == 999;
        censored        = lCensored | rCensored;
        
        meanConvtime    = mean([stackedYl stackedYu], 2);
        
        convtimeRange   = linspace(-72, 120, 5);
        meanBintime     = mean([convtimeRange(1:end-1); ...
             convtimeRange(2:end)], 1);
        
        binIdx          = discretize(meanConvtime, convtimeRange);
    end
    states          = mat2cell(stackedStates, counts);
end

K               = nStates;
[n, d]          = size(stackedPhi);

initDist        = normalize(accumarray(stackedStates, 1));
initPrior       = ones(size(initDist));

transMat        = mkLeftToRight(countTransitions(states, nStates));
transPrior      = triu(ones(nStates), 1);

if strcmpi(cpdType, 'gaussian')
    populationMean      = nanmean(stackedPhi);
    populationSigma     = diag(nanvar(stackedPhi));
    
    mu                  = nan(d, K);
    Sigma               = nan(d, d, K);
    
    for ii=1:nStates
        mask                = stackedStates == ii;
        
        obsPhi              = sum(~isnan(stackedPhi(mask, :)), 1) >= 3;
        
        mu(:, ii)                   = nanmean(stackedPhi(mask, :));
        Sigma(obsPhi, obsPhi, ii)   = nancov(stackedPhi(mask, obsPhi));
        
        mu(~obsPhi, ii)                 = populationMean(~obsPhi);
        Sigma(~obsPhi, ~obsPhi, ii)     = populationSigma(~obsPhi, ~obsPhi);
        
        % check to avoid returning a matrix whose logdet is < eps
        while ~isposdef(Sigma(:, :, ii))
            Sigma(:, :, ii)     = Sigma(:, :, ii) + diag(0.001*ones(d, 1));
        end
    end
    
    cpd             = condGaussCpdCreate(mu, Sigma);
elseif strcmpi(cpdType, 'discrete')
    d                       = size(stackedPhi, 2);

    nObsStates              = nDiscreteValues;
    
    nDiscreteValues         = max(max(stackedPhi(:)), nDiscreteValues);
    tPrior                  = zeros(nStates, nDiscreteValues, d);
    
    for ii=1:d
        columnMax   = max(stackedPhi(stackedPhi(:, ii) <= ...
            nObsStates, ii));
        
        tPrior(:, 1:columnMax, ii)      = 1.1;
        tPrior(:, columnMax+1:end, ii)  = 1;
    end
    
    T           = rand(nStates, nDiscreteValues, size(stackedPhi, 2));
    cpd         = condDiscreteProdCpdCreate(T, 'prior', ...
        struct('alpha', tPrior));
    cpd         = cpd.fitFn(cpd, stackedStates, stackedPhi);
end

initModel.nstates       = K;
initModel.T             = nStates;

initModel.pi            = initDist;
initModel.piPrior       = initPrior;

initModel.A             = transMat;
initModel.transPrior    = transPrior;

initModel.emission      = cpd;

end

function states = getInitStates(stackedData, nStates, cpdType)

if strcmpi(cpdType, 'gaussian')
    mixModel        = mixGaussFit(stackedData, nStates, ...
        'verbose', false, 'maxIter', 50);
    pZ              = mixGaussInferLatent(mixModel, stackedData);
elseif strcmpi(cpdType, 'discrete')
    mixModel        = mixDiscreteFit(stackedData, nStates, ...
        'verbose', false, 'maxIter', 50);
    pZ              = mixDiscreteInferLatent(mixModel, stackedData);
end

[~, states]     = max(pZ, [], 2);

end