function [initModel, states] = getInitModel(data, nStates, cpdType, ...
    varargin)

[   method                                          , ...
    leftToRight                                     , ...
    nDiscreteValues                                 , ...
    modelMissingness                                , ...
    nExpansions                                     , ...
    score                                           , ...
    convtimeLims                                    , ...
    sortModel                                       , ...
    nPerturbGroups                                  , ...
    dx                                              , ...
    yl                                              , ...
    yu                  ] = process_options(varargin, ...
    'method'            , 'mixture'                 , ...
    'leftToRight'       , true                      , ...
    'nDiscreteValues'   , 10                        , ...
    'modelMissingness'  , false                     , ...
    'nExpansions'       , 0                         , ...
    'score'             , {}                        , ...
    'convtimeLims'      , [-72 120]                 , ...
    'sortModel'         , true                      , ...
    'nPerturbGroups'    , 5                         , ...
    'dx'                , true                      , ...
    'yl'                , {}                        , ...
    'yu'                , {}                        );

stackedPhi      = cat(1, data{:});
counts          = cellfun(@(seq)size(seq, 1), data);

if strcmpi(method, 'score') || strcmpi(method, 'score-perturbed')
    assert(~isempty(score), 'Must provide scores for score-based!');
    stackedScore    = cat(1, score{:});
    
    binLimits       = quantile(stackedScore, linspace(0, 1, nStates+1));
    stackedStates   = discretize(stackedScore, binLimits);
    
    if any(isnan(stackedStates))
        index       = arrayfun(@(ii)ones(counts(ii), 1)*ii, ...
            1:length(counts), 'UniformOutput', false);
        index       = cat(1, index{:});
        
        mask            = ~isnan(stackedStates);
        stackedStates   = stackedStates(mask);
        stackedPhi      = stackedPhi(mask, :);
        
        stackedDx       = cat(1, dx{:});
        stackedDx       = stackedDx(mask);
        
        index           = index(mask);
        counts          = histc(index, unique(index));
        
        data            = mat2cell(stackedPhi, counts);
        dx              = mat2cell(stackedDx, counts);
    end
    
    if sortModel
        stackedDx   = cat(1, dx{:});
        [~, order]  = sortHmmModel(nStates, stackedStates, stackedDx, []);

        stackedStates   = colvec(order(stackedStates));
    end
    
    if strcmpi(method, 'score-perturbed')
        groupIdx        = discretize(stackedStates, linspace(0, nStates, ...
            nPerturbGroups+1));
        for gg=rowvec(unique(groupIdx))
            mask        = find(groupIdx == gg);
            stackedStates(mask)     = stackedStates(mask(randperm(...
                length(mask))));
        end
    end
    
    states          = mat2cell(stackedStates, counts);
    initModel       = hmmProfileFitFullyObs(states, data, ...
        cpdType, 'nStates', nStates, 'nDiscreteValues', nDiscreteValues, ...
        'nExpansions', nExpansions, 'modelMissingness', modelMissingness, ...
        'leftToRight', leftToRight);
elseif strcmpi(method, 'convtime')
    assert(~isempty(yl) && ~isempty(yu), ...
        'Must provide convtime for convtime-based!');
    
    stackedYl       = cat(1, yl{:});
    stackedYu       = cat(1, yu{:});

    patCens         = cellfun(@(lower, upper)any(lower == -999) || ...
        any(upper == 999), yl, yu);
    
    meanConvtime    = mean([stackedYl stackedYu], 2);
    
    convtimeBins    = linspace(convtimeLims(1), convtimeLims(2), ...
        nStates+1);
    % flip the states; a high state number should mean a terminal disease
    % stage, but not flipping it would assign a high state number of
    % someone that's actually furthest away from disease!
    stackedStates   = nStates - discretize(meanConvtime, convtimeBins) + 1;
    states          = mat2cell(stackedStates, counts);
    
    if any(patCens)
        initModel   = hmmProfileFitFullyObs(states(~patCens), ...
            data(~patCens), cpdType, 'nStates', nStates, ...
            'nDiscreteValues', nDiscreteValues, ...
            'nExpansions', nExpansions, 'modelMissingness', ...
            modelMissingness);
    
        censIdx     = find(patCens);
        for ii=rowvec(censIdx)
            states{ii}   = colvec(hmmMap(initModel, data{ii}'));
        end
    else
        initModel       = hmmProfileFitFullyObs(states, data, cpdType, ...
            'nStates', nStates, 'nDiscreteValues', nDiscreteValues, ...
            'nExpansions', nExpansions, 'modelMissingness', ...
            modelMissingness, 'leftToRight', leftToRight);
    end
elseif strcmpi(method, 'mixture')
    if strcmpi(cpdType, 'gaussian')
        mixModel        = mixGaussFit(stackedPhi, nStates, ...
            'verbose', false, 'maxIter', 50);
        pZ              = mixGaussInferLatent(mixModel, stackedPhi);
    elseif strcmpi(cpdType, 'discrete')
        mixModel        = mixDiscreteFit(stackedPhi, nStates, ...
            'verbose', false, 'maxIter', 50);
        pZ              = mixDiscreteInferLatent(mixModel, stackedPhi);
    end
    
    [~, stackedStates]  = max(pZ, [], 2);
    
    if sortModel
        stackedDx   = cat(1, dx{:});
        order       = sortHmmModel(nStates, stackedStates, stackedDx, []);

        stackedStates   = order(stackedStates);
    end
    
    states              = mat2cell(stackedStates, counts);
    
    initModel       = hmmProfileFitFullyObs(states, data, ...
        cpdType, 'nStates', nStates, 'nDiscreteValues', nDiscreteValues, ...
        'nExpansions', nExpansions, 'modelMissingness', ...
        modelMissingness, 'leftToRight', leftToRight);
elseif strcmpi(method, 'random')
    stackedStates   = randi(nStates, size(stackedPhi, 1), 1);
    
    states          = mat2cell(stackedStates, counts);
    
    % Get a random placeholder so all the fields, etc. are filled in
    initModel       = hmmProfileFitFullyObs(states, data, ...
        cpdType, 'nStates', nStates, 'nDiscreteValues', nDiscreteValues, ...
        'nExpansions', nExpansions, 'modelMissingness', ...
        modelMissingness, 'leftToRight', leftToRight);
    
    initModel.pi    = normalize(rand(size(initModel.pi)));
    initModel.A     = mkLeftToRight(rand(size(initModel.A)));
    
    T                               = rand(size(initModel.emission.T));
    T(initModel.emission.prior.alpha == 1)  = 0;
    T                               = normalize(T, 2);
    initModel.emission.T            = T;
end

end