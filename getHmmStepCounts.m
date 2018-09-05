function stepCount = getHmmStepCounts(model, startState, endState, varargin)

[   sampleTillReached                                   , ...
    nSamples                                            , ...
    disableBackwardTrans                                , ...
    T                       ] = process_options(varargin, ...
    'sampleTillReached'     , false                     , ...
    'nSamples'              , 100                       , ...
    'disableBackwardTrans'  , true                      , ...
    'T'                     , 25                        );

stepCount       = nan(nSamples, 1);

transMat        = model.A;
if disableBackwardTrans
    transMat    = triu(transMat, 0);
end

K               = size(transMat, 1);

if sampleTillReached
    T           = inf;
end

for s=1:nSamples
    done        = false;
    stateSeq    = startState;
    
    while ~done
        nextState       = randsample(1:K, 1, true, ...
            transMat(stateSeq(end), :));
        stateSeq        = cat(2, stateSeq, nextState);
        
        done            = (nextState == endState) || (length(stateSeq) == T);
    end
    
    if stateSeq(end) == endState
        stepCount(s)    = length(stateSeq)-1;
    else
        stepCount(s)    = T+1;
    end
end

end