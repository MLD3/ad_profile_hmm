function plotTrajectories(stateCountMatrix, stateTransMatrix, subNodes, ...
    colorScaler, varargin)
%% Plot a transition diagram
% stateCountMatrix  -- nNodes x T matrix
%   each value contains the number of instances 

[   T                                                   , ...
    minTrans                                            , ...
    nodeLabels                                          , ...
    censColour                                          , ...
    arrowColour             ] = process_options(varargin, ...
    'T'                     , 10                        , ...
    'minTrans'              , 0.05                      , ...
    'nodeLabels'            , {}                        , ...
    'censColour'            , [0.8903 0.3480 0.1983]    , ...
    'arrowColour'           , [0.8903 0.3480 0.1983]    );

nStates             = size(stateCountMatrix, 1)-1;
nNodes              = length(subNodes);

fullCountMat        = stateCountMatrix(subNodes, 1:T);
fullTransMat        = stateTransMatrix(subNodes, subNodes, 1:T-1);

excludeNodes        = ~ismember(1:nStates+1, subNodes);

excludeCensored     = cumsum(squeeze(sum(stateTransMatrix(...
    excludeNodes, end, :))));

% Remove censored patients that came from stages we don't consider
fullCountMat(subNodes == nStates+1, 2:end)      = ...
    fullCountMat(subNodes == nStates+1, 2:end) - ...
    rowvec(excludeCensored(1:T-1));

fullProbMat         = normalize(fullCountMat, 1);
maxXWidth           = 0.15;
subXWidth           = 0.1;

totalHeight         = 1;
height              = normalize(fullProbMat, 1)*totalHeight;

barStartPoints      = nan(nNodes, T);

if isempty(nodeLabels)
    nodeLabels      = arrayfun(@num2str, subNodes, 'uniformoutput', false);
end

MP                  = get(0, 'MonitorPositions');

fig = figure; hold on;
cbar = colorbar;
fig.Units           = 'pixels';
fig.OuterPosition   = MP(1, :);
drawnow;
ax                  = gca;
ax.Units            = 'inches';
drawnow;
ax.YLim             = [0 totalHeight];
ax.XLim             = [0.7 T+0.3];
% ax.Position         = [0.1 0.6 26.4 13];
% ax.Position         = [0.1 0.6 19.85 10.5];
ax.Position         = [0.58 1.18 18.2 9.80];
ax.YTick            = [];
ax.FontSize         = 30;
ax.YAxis.Visible    = 'on';
ax.TickLength       = [0 0];
ax.XTickLabel       = (ax.XTick-1)*6;
ax.XLabel.String    = 'Time (months)';
ax.YLabel.String    = 'Disease stage';
drawnow;

cmap                = viridis;
hsv                 = rgb2hsv(cmap);
colours             = interp1(linspace(0, 1, size(cmap,1)), hsv, ...
    colorScaler);
colours             = hsv2rgb(colours);

colours             = [censColour; colours];

bars                = cell(nNodes, T);

barAlpha            = 0.8;
edgeAlpha           = 0.80;

% Draw the bars

for t=1:T
    currentY        = 0;
    
    for nn=1:nNodes
        hh          = height(nn, t);
        
        bb          = fill([t-maxXWidth t+maxXWidth t+maxXWidth ...
            t-maxXWidth], [currentY currentY currentY+hh ...
            currentY+hh], colours(nn, :));
        bb.EdgeColor    = bb.FaceColor;
        
        if nn < nNodes
            bb.FaceAlpha    = barAlpha;
        else
            bb.FaceAlpha    = 0.7;
        end
        bb.EdgeAlpha    = 0.1;
        
        bars{nn, t}             = bb;
        barStartPoints(nn, t)   = currentY;
        
        if (t == 1) && (subNodes(nn) ~= (nStates+1))
            text(mean([t-maxXWidth t+maxXWidth]), ...
                currentY + hh/2, nodeLabels{nn}, ...
                'fontsize', 28, 'horizontalalignment', 'center', ...
                'verticalalignment', 'middle');
        end
        
        
        if t < T && subNodes(nn) ~= (nStates+1)
            transProb       = normalize(fullTransMat(:, :, t), 2);
            
            selfProb        = transProb(nn, nn);
            censProb        = transProb(nn, subNodes == nStates+1);
            
            xx          = [t+maxXWidth-subXWidth t+maxXWidth...
                t+maxXWidth t+maxXWidth-subXWidth];
            
            censBar         = fill(xx, ...
                [currentY+selfProb*hh currentY+selfProb*hh ...
                currentY+(selfProb+censProb)*hh ...
                currentY+(selfProb+censProb)*hh], [1 1 1]);
            
            censBar.EdgeAlpha   = 0;
            
            censBar         = fill(xx, ...
                [currentY+selfProb*hh currentY+selfProb*hh ...
                currentY+(selfProb+censProb)*hh ...
                currentY+(selfProb+censProb)*hh], censColour);
            
            censBar.EdgeAlpha   = 0;
            censBar.FaceAlpha   = 0.8;
        end
        
        currentY        = currentY + hh;
    end
end

edges           = cell(nNodes, T-1);
maxEdgeWidth    = 0.03;
maxArrWidth     = 72*0.2;
maxArrLength    = 72*0.2;

xScaler         = ax.Position(3)/(ax.XLim(2)-ax.XLim(1));
yScaler         = ax.Position(4)/(ax.YLim(2)-ax.YLim(1));

for t=1:T-1
    transMat        = normalize(fullTransMat(:, :, t), 2);
    
    for src=1:nNodes
        maxH            = height(src, t);
        sourceY         = barStartPoints(src, t) + ...
            (transMat(src, src) + transMat(src, ...
            subNodes == nStates+1))*maxH;
        
        for dest=1:nNodes
            if (dest ~= src) && (transMat(src, dest) >= minTrans) && ...
                    (dest ~= find(subNodes == nStates+1))
                arrowH      = transMat(src, dest)*maxH;
                
                destY       = barStartPoints(dest, t+1) + ...
                    height(dest, t+1)/2 - (arrowH/2);
                
                dY          = (destY-sourceY);
                dX          = (t+1-maxXWidth) - (t+maxXWidth);
                
                angle       = atand((dY*yScaler)/(dX*xScaler));
                hh          = arrowH; %/cosd(angle);
                
                ee      = fill([t+maxXWidth t+0.5 ...
                    t+0.5 t+maxXWidth], [sourceY+hh ...
                    (sourceY+dY/2+hh) (sourceY+dY/2) sourceY], ...
                    bars{src, t}.FaceColor);
                ee.FaceAlpha        = edgeAlpha;
                ee.EdgeColor        = ee.FaceColor;
                ee.EdgeAlpha        = 0;
                
                ee      = fill([t+0.5 t+1-maxXWidth ...
                    t+1-maxXWidth t+0.5], [(sourceY+dY/2+hh) ...
                    destY+hh destY (sourceY+dY/2)], ...
                    bars{dest, t+1}.FaceColor);
                ee.FaceAlpha        = edgeAlpha;
                ee.EdgeColor        = ee.FaceColor;
                ee.EdgeAlpha        = 0;
                
                edges{src, t}       = ee;
                
                ar                  = annotation('arrow');
                ar.Units            = 'inches';
                ar.X        = [ax.Position(1)+(t-ax.XLim(1)+...
                    maxXWidth)*xScaler ...
                    ax.Position(1)+(t-ax.XLim(1)+0.5)*xScaler];
                ar.Y        = [ax.Position(2) + ...
                    (sourceY+hh/2-ax.YLim(1))*yScaler ...
                    ax.Position(2) + ...
                    (sourceY-ax.YLim(1)+(dY/2)+hh/2)*yScaler];
                ar.LineStyle        = 'none';
                ar.HeadWidth        = maxArrWidth;
                ar.HeadLength       = maxArrLength;
                ar.Color            = arrowColour;
                
                sourceY             = sourceY+hh;
            end
        end
    end
end

colormap(cmap);
cbar.Ticks = [0 0.2 0.4 0.6 0.8 1.0];
cbar.TickLabels{1} = 'AD';
cbar.TickLabels{end} = 'MCI';

end