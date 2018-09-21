function [tform2, warpErrors] = Warp2BestPair(hybe1,hybe2,varargin)
% Compute translation/rotation warp that best aligns the points in hybe1
% and hybe2 by maximizing the alignment of the two points that show the 
% most mutually consistent x,y translation.  
% Copyright Presidents and Fellows of Harvard College, 2016.


% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'maxD', 'nonnegative', 2};
defaults(end+1,:) = {'useCorrAlign', 'boolean', true};
defaults(end+1,:) = {'fighandle', 'handle', []};
defaults(end+1,:) = {'imageSize', 'array', [256 256]};
defaults(end+1,:) = {'showPlots', 'boolean', true};
defaults(end+1,:) = {'verbose', 'boolean', true};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabSTORM:invalidArguments', 'two nx2 vectors of points are required');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);
% parameters = ParseVariableArguments([], defaults);

%--------------------------------------------------------------------------
%% Main Function
%--------------------------------------------------------------------------
[matched1, matched2,parameters] = MatchFeducials(hybe1,hybe2,'parameters',parameters);

if length(matched1) < 2
    error('found fewer than 2 feducials, cannot compute warp');
else

    % Warp that maximizes MSD for best pair of beads
    % Compare the shift vectors, identify the pair that have the most similar 
    shifts = [hybe1(matched1,1) - hybe2(matched2,1), hybe1(matched1,2) - hybe2(matched2,2)];
    [idx,dist] = knnsearch(shifts,shifts,'K',2);
    [~,shiftIdx]  = min(dist(:,2));
    bestpair = [matched1(shiftIdx),matched2(shiftIdx);
            matched1(idx(shiftIdx,2)),matched2(idx(shiftIdx,2))];
    tform2 = WarpPoints(hybe1(bestpair(:,1),:),hybe2(bestpair(:,2),:),'translation rotation');
    [xw,yw] = tforminv(tform2,hybe2(:,1),hybe2(:,2));

    if parameters.showPlots
        if isempty(parameters.fighandle);
            figure;  
        else
            figure(parameters.fighandle); 
        end
        plot(xw,yw,'k.'); hold on;
        plot(hybe1(:,1),hybe1(:,2),'bo');
        title('aligned by the pair of beads that best match');
         pause(.2);
         set(gcf,'color','w');
    end

    [~,dist_PairWarp] = knnsearch(hybe1,[xw,yw]);
    tform2.tdata.Tinv(1,1) = 1;
    dist_PairWarp(dist_PairWarp > parameters.maxD) = NaN;
    warpErrors = sort(dist_PairWarp);

end
