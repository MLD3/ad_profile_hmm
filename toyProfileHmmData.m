T               = 5;
K               = T*2 + 2;

stateKey        = {'BEGIN', 'MATCH', 'DELETE', 'END'};
stateType       = [1 repmat([2 3], 1, T) 4];

matchStates             = find(stateType == 2);
deleteStates            = find(stateType == 3);

A               = zeros(K);
A(1, [2 3])     = 1;
A(end, end)     = 0;

for ii=1:T
    matchIdx    = matchStates(ii);
    % Can transition from a match state to the next match/delete state 
    A(matchIdx, min(matchIdx+[2 3], K))     = 1;
    % allow for self-transitions with half the likelihood
    A(matchIdx, matchIdx)   = 0.5;
    
    delIdx      = deleteStates(ii);
    % Can transition from a delete state to the next match/delete state 
    A(delIdx, min(delIdx+[1 2], K))         = 1;
end

A           = normalize(A, 2);

d           = T;
mu          = nan(d, K);
Sigma       = nan(d, d, K)*eps;

mu(:, matchStates)      = ones(d, T)*0.01;
diagIdx                 = sub2ind([d K], 1:d, matchStates);
mu(diagIdx)             = 0.96;

Sigma(:, :, matchStates)    = repmat(eye(d), 1, 1, T);

emissionType    = 'multinomial';

nSamples        = 100;
observations    = cell(nSamples, 1);
states          = cell(nSamples, 1);
likelihood      = cell(nSamples, 1);

for ii=1:nSamples
    currentState        = 1;
    
    obsSeq              = [];
    stateSeq            = [currentState];
    likelihoodSeq       = [1];
    
    while currentState ~= K
        nextState       = randsample(1:K, 1, true, A(currentState, :));
        likelihoodSeq   = cat(1, likelihoodSeq, ...
            likelihoodSeq(end)*A(currentState, nextState));
        
        if ismember(nextState, matchStates)
            if strcmpi(emissionType, 'multinomial')
                obsSeq      = cat(1, obsSeq, ...
                    randsample(1:d, 1, true, mu(:, nextState)));
                likelihoodSeq(end)      = likelihoodSeq(end)*...
                    mu(obsSeq(end), nextState);
            elseif strcmpi(emissionType, 'gaussian')
                obsSeq      = cat(1, obsSeq, ...
                    mvnrnd(mu(:, nextState), Sigma(:, :, nextState), 1));
                likelihoodSeq(end)      = likelihoodSeq(end)*...
                    mvnpdf(obsSeq(end, :), rowvec(mu(:, nextState)), ...
                    Sigma(:, :, nextState));
            end
        end
        
        stateSeq        = cat(2, stateSeq, nextState);
        currentState    = nextState;
    end
    
    observations{ii}    = obsSeq;
    states{ii}          = stateSeq;
    likelihood{ii}      = likelihoodSeq;
end

obsLen              = cellfun(@(seq)size(seq, 1), observations);
nontrivial          = obsLen > 0;

observations        = observations(nontrivial);
states              = states(nontrivial);
likelihood          = likelihood(nontrivial);

initModel.nstates       = K;
initModel.T             = T;
initModel.A             = A;
initModel.emission      = struct('mu', mu, 'Sigma', Sigma, 'type', ...
    'multinomial');
initModel.stateType     = stateKey(stateType);

save('profile-hmm-toy-data.mat');
