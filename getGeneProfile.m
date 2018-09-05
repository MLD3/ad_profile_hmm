function geneProfile = getGeneProfile(data)

if istable(data)
    apgen           = data{:, {'CAT_APGEN1_2_0', ...
        'CAT_APGEN1_3_0', 'CAT_APGEN1_4_0', 'CAT_APGEN2_2_0', ...
        'CAT_APGEN2_3_0', 'CAT_APGEN2_4_0'}};
elseif ismatrix(data)
    apgen           = data;
end
    
[apgen1, a1]    = find(apgen(:, 1:3)');
[apgen2, a2]    = find(apgen(:, 4:6)');

allele          = nan(size(apgen, 1), 2);
allele(a1, 1)   = apgen1+1;
allele(a2, 2)   = apgen2+1;

geneProfile     = nan(size(allele, 1), 1);

for ii=1:length(geneProfile)
    if isequal(allele(ii, :), [2 2])
        geneProfile(ii)     = 1;
    elseif isequal(allele(ii, :), [3 3])
        geneProfile(ii)     = 2;
    elseif isequal(allele(ii, :), [4 4])
        geneProfile(ii)     = 3;
    elseif isequal(allele(ii, :), [2 3])
        geneProfile(ii)     = 4;
    elseif isequal(allele(ii, :), [3 4])
        geneProfile(ii)     = 5;
    elseif isequal(allele(ii, :), [2 4])
        geneProfile(ii)     = 6;
    end
end