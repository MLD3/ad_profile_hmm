function [Y, allSamples] = getSampleQuantiles(cpd, discretePhi, rawPhi, ...
    stateIdx, stackedStates, varargin)

[   nBootSamples                                        , ...
    quantileVals                                        , ...
    sampleObs               ] = process_options(varargin, ...
    'nBootSamples'          , 10000                     , ...
    'quantileVals'          , [0.025 0.25 0.5 0.75 0.975], ...
    'sampleObs'             , false                     );

T               = cpd.T;
nObsStates      = size(T, 2);
nStates         = length(stateIdx);

d               = size(discretePhi, 2);
Y               = zeros(length(quantileVals), nStates, d);
allSamples      = zeros(nBootSamples, nStates, d);

for ii=1:d
    rawValues       = rawPhi(:, ii);
    discValues      = discretePhi(:, ii);
    
    isObserved      = ~isnan(discValues);
    
    for k=1:nStates
        if sampleObs
            weights         = T(stateIdx(k), discValues(isObserved), ii);
            samples         = randsample(rawValues(isObserved), ...
                nBootSamples, true, weights);
            
            % mnSamples       = mnrnd(1, T(stateIdx(k), :, ii), nBootSamples);
            % [discreteSamples, ~]    = find(mnSamples');
            % samples                 = nan(size(discreteSamples));
            %
            % for o=1:nObsStates
                % obsMask             = discreteSamples == o;
                % if sum(discValues == o) > 0
                %   samples(obsMask)    = randsample(...
                %       rawValues(discValues == o), sum(obsMask), true);
                % else
                %   samples(obsMask)    = nan;
                % end
            % end
            
            allSamples(:, k, ii)    = samples;
            Y(:, k, ii)             = quantile(samples, quantileVals);
        else
            stateMask       = stackedStates == k;
            if nBootSamples > 1
                samples     = randsample(rawValues(stateMask), ...
                    nBootSamples, true);
            else
                samples     = rawValues(stateMask);
            end
            
            Y(:, k, ii)     = quantile(samples, quantileVals);
        end
    end
end