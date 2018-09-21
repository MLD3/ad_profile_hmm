function [phi, phiSum, summaryNames]  = getClinicalSummary(data, varargin)

[   scale                                   , ...
    imputeNan                               , ...
    phiNames                                , ...
    counts                                  , ...
    age         ] = process_options(varargin, ...
    'scale'     ,   true                    , ...
    'imputeNan' ,   true                    , ...
    'phiNames'  ,   {}                      , ...
    'counts'    ,   []                      , ...
    'age'       ,   true                    );

if ~istable(data)
    assert(~isempty(phiNames), ...
        'Must provide features names to retrieve summary');
end

% ADAS-Cog, TRAILS-B, FAQTOTAL, RAVLT, MMSE, CDR

summaryNames    = {'ADAS', 'TRAILS', 'FAQ', 'RAVLT', 'MMSE', 'CDR'};

phiCodes        = {{'Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6', 'Q7', 'Q8', ...
    'Q9', 'Q10', 'Q11', 'Q12', 'Q13'}, {'TRABSCOR'}, ...
    {'FAQFINAN', 'FAQFORM', 'FAQSHOP', 'FAQGAME', 'FAQBEVG', ...
    'FAQMEAL', 'FAQEVENT', 'FAQTV', 'FAQREM', 'FAQTRAVL'}, ...
    {'AVTOT1', 'AVTOT2', 'AVTOT3', 'AVTOT4', 'AVTOT5', 'AVTOT6'}, ...
    {'MMORIENTATION', 'MMRECALL', 'MMATTENTION', 'MMDLRECALL', ...
    'MMLANGUAGE', 'MMCOMMAND', 'MMREPEAT', 'MMREAD', 'MMWRITE', 'MMDRAW'}, ...
    {'CDMEMORY', 'CDORIENT', 'CDJUDGE', 'CDCOMMUN', 'CDHOME', ...
    'CDCARE', 'CDGLOBAL'}};

phiSum      = zeros(size(data, 1), length(phiCodes));

for testIdx=1:length(phiCodes)
    test        = phiCodes{testIdx};
    for subscore = test
        if istable(data)
            col     = data{:, subscore};
        elseif ismatrix(data)
            col     = data(:, strcmpi(phiNames, subscore));
        end
        
        if imputeNan
            col(isnan(col))     = nanmean(col);
        end
        
        phiSum(:, testIdx)   = phiSum(:, testIdx) + col;
    end
end

if age
    summaryNames    = cat(2, summaryNames, 'AGE');
    if istable(data)
        phiSum      = [phiSum data.AGE];
    elseif ismatrix(data)
        phiSum      = [phiSum data(:, strcmpi(phiNames, 'AGE'))];
    end
end

if scale
    phiSum      = zScale(phiSum);
end

if ~isempty(counts)
    phi         = mat2cell(phiSum, counts);
else
    phi         = phiSum;
end

end