function sp = dealSubplots(nRows, nCols)
%function sp = dealSubplots(nRows, nCols)
%
% Generates a large number of subplots within a figure for plotting many
% pannels easily. % Does some basic formatting for tidyness (e.g. smaller fontsize, no
% background color).
%
% INPUTS:
%   nRows: number of rows
%   nColumns: number of columns
%
% OUTPUTS:
%   sp: matrix (nRows x nColumns) containing axes handles
%
% Stephen Town - 2015

% Preassign
nPlots = nRows*nCols;
sp = nan(nPlots,1);

% Create axes
for i = 1 : nPlots
    sp(i) = subplot(nRows, nCols, i);
end

% Format axes
set(sp,'nextPlot','add',...
    'FontSize',8,...
    'color','none')
