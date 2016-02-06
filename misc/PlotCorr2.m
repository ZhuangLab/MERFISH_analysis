function pearsonCorr = PlotCorr2(x,y,varargin)
% PlotCorr2(x,y)
% creates log-log plot of correlation betweent x and y
% stats = PlotCorr(x,y) returns stats.log10rho, stast.log10pvalue
%  stats.rho and stats.pvalue for the correlation in addition to the plot
%
%  non-zero points are removed from log-log correlation and correlation
%  plot, but not from the linear correlations computed.  
% Second version of PlotCorr
% Copyright Presidents and Fellows of Harvard College, 2016.

% -------------------------------------------------------------------------
% Default variables
% -------------------------------------------------------------------------
defaults = cell(0,3);
defaults(end+1,:) = {'MarkerSize', 'positive', 10};
defaults(end+1,:) = {'FontSize', 'positive', 6};
defaults(end+1,:) = {'colorMap', 'colormap', 'jet'};
defaults(end+1,:) = {'nameBuffer', 'positive', .1};
defaults(end+1,:) = {'figHandle', 'handle', []};
defaults(end+1,:) = {'axesHandle', 'handle', []};
defaults(end+1,:) = {'pointNames', 'cell', {}};
defaults(end+1,:) = {'plotFunction', 'function', @loglog};
defaults(end+1,:) = {'includeLog10', 'boolean', true};
defaults(end+1,:) = {'includeLin', 'boolean', true};

% -------------------------------------------------------------------------
% Parse necessary input
% -------------------------------------------------------------------------
if nargin < 2
    error('matlabSTORM:invalidArguments', 'requires x,y');
end

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
parameters = ParseVariableArguments(varargin, defaults, mfilename);

% -------------------------------------------------------------------------
% Parse variable input
% -------------------------------------------------------------------------
non0 = x>0 & y>0;
nonNaN = ~(isnan(x) | isnan(y) | isinf(x) | isinf(y) );

if ~isempty(parameters.figHandle)
    set(0, 'CurrentFigure', parameters.figHandle);
end
if ~isempty(parameters.axesHandle)
    axes(parameters.axesHandle)
end

if isempty(x(non0 & nonNaN));
   warning('no nonzero data to plot'); 
   pearsonCorr.log10rho = [];
   pearsonCorr.log10pvalue = [];
   pearsonCorr.rho = [];
   pearsonCorr.pvalue = []; 
else
    if ~isempty(parameters.pointNames);
        parameters.plotFunction(x,y,'k.','MarkerSize',parameters.MarkerSize); hold on;
        text(x+parameters.nameBuffer*x,y,...
            parameters.pointNames,'FontSize',parameters.FontSize); 
    else
        parameters.plotFunction(x,y,'k.','MarkerSize',parameters.MarkerSize); hold on;
    end
    x= reshape(x, [numel(x) 1]);
    y = reshape(y, [numel(y) 1]);
    
    % Calculate correlation coefficient
    [c1,p1] = corr(x(nonNaN),y(nonNaN));
    pearsonCorr.rho = c1;
    pearsonCorr.pvalue = p1; 
    titleString = '';
    % Calculate correlation coefficient for log10 data
    if parameters.includeLog10
        [c0,p0] = corr(log10(x(non0 & nonNaN)),log10(y(non0 & nonNaN)));

        pearsonCorr.log10rho = c0;
        pearsonCorr.log10pvalue = p0;
        titleString = ['\rho_{log10} = ',num2str(c0,2),' (p=',num2str(p0,2),')  '];
    end
    if parameters.includeLin
        titleString = [titleString '\rho = ',num2str(c1,2),' (p=',num2str(p1,2),')'];
    end
    title(titleString);
end