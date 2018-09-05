function kld = getMeanKld(freqMat)

probHist        = normalize(freqMat+1, 1);

allPairs        = nchoosek(1:size(probHist, 2), 2);

kld             = 0;

for ii=1:size(allPairs, 1)
    pair        = allPairs(ii, :);
    
    h1          = probHist(:, pair(1));
    h2          = probHist(:, pair(2));
    
    kld         = kld + ...
        mean([sum(h1 .* log((h1)./(h2))) sum(h2 .* log(h2./h1))]);
end

end