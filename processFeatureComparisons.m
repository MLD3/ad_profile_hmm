counts              = allData.counts;

countGroups         = unique(counts);
groupPop            = histc(counts, countGroups);

validGroups         = find(groupPop >= 10);
nValidGroups        = length(validGroups);
pValMatrix          = nan(nHmmExp, nHmmExp, nValidGroups);

flatData            = reshape(loglikZ, size(loglikZ, 1), []);

for ii=1:nValidGroups
    groupMask       = counts == countGroups(validGroups(ii));
    
    pValMatrix(:, :, ii)    = pairwiseHypothesisTest(...
        flatData(groupMask, :), 'wilcoxon');
end

meanPval        = mean(pValMatrix, 3);


