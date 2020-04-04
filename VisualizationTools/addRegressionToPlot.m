function [mdl, h_out] = addRegressionToPlot(ax, opt)
% function [mdl, p] = addRegressionToPlot(ax, opt)
%
% Fits a linear regression to data points within a defined axis. 
%
% INPUTS (required):
%   ax: handles for  the axis containing data points 
%
% INPUTS (optional): 
%   opt: visualization options including whether to replot the regression function
%   confidence intervals, legend and if data points should be replotted (an
%   annoying feature of matlab)
%
%   e.g opt = struct('plotData',true,'plotFit',true,'plotCI',true,...
%                  'display',true,'legend',true);
%
% OUTPUTS:
%   mdl: linear model object
%   h_out: handles of objects plotted
%
% Stephen Town - 2015

% (To do in future - enable handles for data points to be enterred into
% function directly. Likewise enable direct input of x,y values

% Default options
if nargin == 1
    opt = struct('plotData',true,'plotFit',true,'plotCI',true,...
                 'display',true,'legend',true);
end

% Get axes properties
if isstruct(ax)
    t = ax.Title.String;    
    xstr = ax.XLabel.String;
    ystr = ax.YLabel.String;
else
    th  = get(ax,'Title');  t    = get(th,'String');
    xsh = get(ax,'XLabel'); xstr = get(xsh,'String');
    ysh = get(ax,'YLabel'); ystr = get(ysh,'String');
end

% Get scatter data
s = findobj(ax,'type','scatter');
NS = numel(s);
[xs, ys, cs] = deal( cell(1,NS));

for i = 1 : NS
    xs{i} = get(s(i),'xdata');
    ys{i} = get(s(i),'ydata');
    cs{i} = get(s(i),'MarkerEdgeColor');
end

% Get line data
h = findobj(ax,'type','line');
NL = numel(h);
[xh, yh, ch] = deal( cell(1,NL));

for i = 1 : NL
    xh{i} = get(h(i),'xdata');
    yh{i} = get(h(i),'ydata');
    ch{i} = get(h(i),'color');
end

% Bring together data
x = [xs xh];
y = [ys yh];
c = [cs ch];

x(cellfun(@isempty, x)) = [];
y(cellfun(@isempty, y)) = [];
c(cellfun(@isempty, c)) = [];
n = numel(x);

% Manage axes
set(ax,'nextPlot','add')

% For each object
for i = 1 : n

    % Fit linear regression
    mdl = fitlm(x{i},y{i});
    % display(mdl)
    
    % Plot model
    h_out = plot(mdl,'parent',ax);
    
    % Parse graphic object names
    for j = 1 : numel(h_out)-1        
        eval(sprintf('pH.%s = p(j);', strrep(h_out(j).DisplayName,' ','')))
    end
    
    % Set plot
    set(h_out(3:4),'LineStyle','--')
        
    % Delete any added model plots that aren't required
    axes(ax)
    if ~opt.plotData, delete(pH.Data); end
    if ~opt.plotFit, delete(pH.Fit); end
    if ~opt.plotCI, delete(h_out(3:4)); end % This will only delete one of the two bounds
    if ~opt.legend, legend off; end
end

% Recover axes labels and titles
xlabel(ax,xstr)
ylabel(ax,ystr)
title(ax,t)
