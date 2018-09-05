function model = getConvtimeModel(data, nStates, cpdType, yl, yu, varargin)
%% Build a simple baseline model that assigns states based on convtime
% Do not use temporal information, this is a baseline model

[   convtimeLims                                        , ...
    nDiscreteValues         ] = process_options(varargin, ...
    'convtimeLims'          , [-72 120]                 , ...
    'nDiscreteValues'       , 10                        );

assert(~isempty(yl) && ~isempty(yu), ...
    'Must provide convtime for convtime-based!');
    
stackedYl       = cat(1, yl{:});
stackedYu       = cat(1, yu{:});
stackedPhi      = cat(1, data{:});

censored        = (stackedYl == -999) | (stackedYu == 999);
    
meanConvtime    = mean([stackedYl stackedYu], 2);
    
convtimeBins    = linspace(convtimeLims(1), convtimeLims(2), ...
    nStates+1);
% flip the states; a high state number should mean a terminal disease
% stage, but not flipping it would assign a high state number of
% someone that's actually furthest away from disease!
stackedStates   = nStates - discretize(meanConvtime, convtimeBins) + 1;
    
if any(censored)
    model       = convtimeFitFullyObs(stackedStates(~censored), ...
        stackedPhi(~censored, :), cpdType, nStates, nDiscreteValues);
    
    if strcmpi(cpdType, 'discrete')
        ss      = mixDiscreteInferLatent(model, stackedPhi);
    elseif strcmpi(cpdType, 'gaussian')
        ss      = mixGaussInferLatent(model, stackedPhi);
    else
        error('Unknown CPD type!');
    end

    [~, stackedStates(censored)]    = max(ss(censored), [], 2);
end

model       = convtimeFitFullyObs(stackedStates, stackedPhi, cpdType, ...
    nStates, nDiscreteValues);

model.bins  = convtimeBins;

end

function model = convtimeFitFullyObs(states, phi, cpdType, nStates, ...
    nDiscreteValues)

if strcmpi(cpdType, 'discrete')
    d               = size(phi, 2);
    
    nObsStates      = nDiscreteValues;
    
    nDiscreteValues         = max(max(phi(:)), nDiscreteValues);
    tPrior                  = zeros(nStates, nDiscreteValues, d);
    
    for ii=1:d
        columnMax   = max(phi(phi(:, ii) <= nObsStates, ii));
        
        tPrior(:, 1:columnMax, ii)      = 2;
        tPrior(:, columnMax+1:end, ii)  = 1;
    end
    
    T           = rand(nStates, nDiscreteValues, size(phi, 2));
    cpd         = condDiscreteProdCpdCreate(T, 'prior', ...
        struct('alpha', tPrior));
    cpd         = cpd.fitFn(cpd, states, phi);
else
    error('Unsupported CPD type!');
end

model       = struct('cpd', cpd, 'nstates', nStates, 'nmix', nStates, ...
    'mixWeight', normalize(accumarray(states, 1, [nStates 1])));

end
