function Y = getSampleExpectations(cpd, qLims, varargin)

[   stateIdx                ] = process_options(varargin, ...
    'stateIdx'              , 1:size(cpd.T, 1)          );

T               = cpd.T;
d               = size(T, 3);
nStates         = length(stateIdx);

Y               = zeros(nStates, d);

for ii=1:d
    binLims         = qLims{ii};
    
    if ~isempty(binLims)
        probMat         = T(:, :, ii);
        
        nObsStates      = length(binLims)-1;
        
        lowerBound      = probMat(:, 1:nObsStates)*colvec(binLims(1:end-1));
        upperBound      = probMat(:, 1:nObsStates)*colvec(binLims(2:end));
        
        Y(:, ii)        = mean([lowerBound upperBound], 2);
    end
end

end