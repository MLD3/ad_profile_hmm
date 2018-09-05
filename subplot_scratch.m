data        = rand(2, 10, 4);

nPlots      = size(data, 3);
fig         = figure;
subs        = cell(nPlots, 1);

for ii=1:nPlots
    hSub        = subplot(nPlots, 1, ii);
    subs{ii}    = hSub;
    
    errorbar(data(1, :, ii), data(2, :, ii), 'o');
    grid on; box on;
    
    if ii < nPlots
        hSub.XTickLabel     = [];
    end
    
    hSub.Position(4)    = hSub.Position(4) + 0.03;
    
    if mod(ii, 2) == 0
        hSub.YAxisLocation  = 'right';
    end
end

linkaxes(cat(2, subs{:}), 'x');

hSub.XLim       = [0.7 10.3];

