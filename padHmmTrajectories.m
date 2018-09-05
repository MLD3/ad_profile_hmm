function [trajPad, dxPad, ylPad, yuPad, agePad, followPad, scorePad]    = ...
    padHmmTrajectories(trajectories, dx, yl, yu, age, followup, score, res)
%% Pad the HMM trajectories to account for uneven sampling rates
% trajectories          -- each cell is ni x d
% followup, yl, yu      -- each cell is ni x 1
% res                   -- the desired minimum time between visits

counts          = cellfun(@length, yl);
[~, d]          = size(trajectories{1});

stackedFollowup     = cat(1, followup{:});
stackedFollowup     = round(stackedFollowup/res)*res;
roundFollowup       = mat2cell(stackedFollowup, counts);

trajPad         = trajectories;
dxPad           = dx;
ylPad           = yl;
yuPad           = yu;
agePad          = age;
scorePad        = score;
followPad       = followup;

for ii=1:length(trajectories)
    followTimes     = roundFollowup{ii};
    followDiffs     = abs(diff(followTimes));
    
    if ~all(followDiffs == res)
        phi             = trajectories{ii};
        ylow            = yl{ii};
        yupper          = yu{ii};
        follow          = followup{ii};
        labels          = dx{ii};
        scr             = score{ii};
        ag              = age{ii};
        
        newPhi          = phi(1, :);
        newYl           = ylow(1);
        newYu           = yupper(1);
        newAge          = ag(1);
        newFollow       = follow(1);
        newLabels       = labels(1);
        newScores       = scr(1);
        
        for jj=1:length(followDiffs)
            if followDiffs(jj) ~= res
                numAdd      = (followDiffs(jj)/res)-1;
                
                newPhi      = cat(1, newPhi, nan(numAdd, d));
                
                if ylow(jj) ~= ylow(jj+1)
                    newYl       = cat(1, newYl, ...
                        colvec(linspace(ylow(jj)-res, ylow(jj+1)+res, ...
                        numAdd)));
                else
                    % if ylow(jj) == ylow(jj+1) == -999
                    newYl       = cat(1, newYl, ones(numAdd, 1)*ylow(jj));
                end
                
                if yupper(jj) ~= yupper(jj+1)
                    newYu       = cat(1, newYu, ...
                        colvec(linspace(yupper(jj)-res, ...
                        yupper(jj+1)+res, numAdd)));
                else
                    % if yupper(jj) == yupper(jj+1) == 999
                    newYu       = cat(1, newYu, ones(numAdd, 1)*yupper(jj));
                end
                
                if labels(jj) ~= labels(jj+1)
                    newLabels   = cat(1, newLabels, nan(numAdd, 1));
                else
                    % if labels(jj) == labels(jj+1), only then interpolate
                    newLabels   = cat(1, newLabels, ...
                        ones(numAdd, 1)*labels(jj));
                end
                
                newFollow   = cat(1, newFollow, ...
                    colvec(linspace(follow(jj)-res, follow(jj+1)+res, ...
                    numAdd)));
                
                newAge      = cat(1, newAge, ...
                    colvec(linspace(ag(jj)-res, ag(jj+1)+res, ...
                    numAdd)));
                
                if follow(jj) == follow(jj+1)
                    newScores   = cat(1, newScores, ones(numAdd, 1)*scr(jj));
                else
                    newScores   = cat(1, newScores, ...
                        colvec(interp1(follow(jj:jj+1), scr(jj:jj+1), ...
                        colvec(linspace(follow(jj)-res, follow(jj+1)+res, ...
                        numAdd)))));
                end
            end
            
            newPhi          = cat(1, newPhi, phi(jj+1, :));
            newYl           = cat(1, newYl, ylow(jj+1));
            newYu           = cat(1, newYu, yupper(jj+1));
            newLabels       = cat(1, newLabels, labels(jj+1));
            newFollow       = cat(1, newFollow, follow(jj+1));
            newAge          = cat(1, newAge, ag(jj+1));
            newScores       = cat(1, newScores, scr(jj+1, :));
        end
        
        trajPad{ii}         = newPhi;
        ylPad{ii}           = newYl;
        yuPad{ii}           = newYu;
        dxPad{ii}           = newLabels;
        agePad{ii}          = newAge;
        scorePad{ii}        = newScores;
        followPad{ii}       = newFollow;
    end
end

matSizeFunc         = @(seq)size(seq, 1);
assert(isequal(cellfun(matSizeFunc, trajPad), cellfun(@length, ylPad), ...
    cellfun(@length, yuPad), cellfun(@length, dxPad), ...
    cellfun(@length, dxPad), cellfun(@length, scorePad), ...
    cellfun(@length, followPad), cellfun(@length, agePad)), ...
    'Length mismatch while padding!');

end