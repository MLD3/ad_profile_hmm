K               = 6;
O               = 10;
d               = 3;

initDist        = normalize(rand(K, 1));
A               = normalize(rand(K), 2);
B               = normalize(rand(K, O, d), 2);

nSamples        = 100;
observations    = cell(nSamples, 1);
states          = cell(nSamples, 1);
likelihood      = cell(nSamples, 1);

for ii=1:nSamples
    currentState        = randsample(K, 1, true, initDist);
    stateSeq            = [currentState];
    
    sample              = mnrnd(1, squeeze(B(currentState, :, :))');
    [obs, ~]            = find(sample');
    obsSeq              = rowvec(obs);
    likelihoodSeq       = [1];
    
    seqLen              = randsample(6, 1);
    
    while length(stateSeq) < seqLen
        nextState       = randsample(1:K, 1, true, A(currentState, :));
        likelihoodSeq   = cat(1, likelihoodSeq, ...
            likelihoodSeq(end)*A(currentState, nextState));
        
        sample          = mnrnd(1, squeeze(B(nextState, :, :))');
        [obs, ~]        = find(sample');
        obsSeq          = cat(1, obsSeq, rowvec(obs));
        likelihoodSeq(end)      = likelihoodSeq(end)*...
            B(nextState, obsSeq(end));
        
        stateSeq        = cat(2, stateSeq, nextState);
        currentState    = nextState;
    end
    
    observations{ii}    = obsSeq;
    states{ii}          = stateSeq;
    likelihood{ii}      = likelihoodSeq;
end

data        = cellfun(@transpose, observations, 'UniformOutput', false);

initModel       = struct('pi', normalize(rand(K, 1)), ...
    'A', normalize(rand(K), 2), 'emission', condDiscreteProdCpdCreate(...
    normalize(rand(K, O, d), 2), 'prior', struct('alpha', 1)), ...
    'nstates', K, 'piPrior', zeros(1, K), 'transPrior', zeros(K));

[model, loglik]         = hmmProfileFitEm(data, 'condDiscreteProd', ...
    'model', initModel, 'EMargs', {'verbose', true});

[modelPmtk, loglikPmtk]     = hmmFitEm(data, K, 'discrete', ...
    'pi0', initModel.pi, 'trans0', initModel.A, 'emission0', ...
    initModel.emission, 'piPrior', initModel.piPrior, ...
    'transPrior', initModel.transPrior, 'emissionPrior', ...
    initModel.emission.prior, 'verbose', true);

logp            = cellfun(@(seq)hmmLogprob(model, seq), data);
logpPmtk        = cellfun(@(seq)hmmLogprob(modelPmtk, seq), data);

X       = data{randsample(length(data), 1)};
B       = exp(mkSoftEvidence(model.emission, X)); 
logB    = mkSoftEvidence(model.emission, X);

[gamma, alpha, beta, logp]  = hmmProfileFwdBack(model.pi, model.A, B);
[gammaPmtk, logpPmtk, alphaPmtk, betaPmtk]  = hmmInferNodes(modelPmtk, ...
    X);
