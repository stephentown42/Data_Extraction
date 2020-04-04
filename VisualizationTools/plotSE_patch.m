function [line_h, patch_h] = plotSE_patch(x,y,s,ax,color)
% [line_h, patch_h] = plotSE_patch(x,y,s,ax,color)
%
% Can be used to calculate mean and standard error if ischar(s) or for
% data on which standard error is already calculated (s = double input)
%
% If calculating standard error, y is matrix with rows as trials and columns
% as time samples.
%
% INPUTS (required):
%   x: vector of values along the x-axis
%   y: either vector of mean values on y axis, or else matrix of individual
%   samples (rows), where number of columns is equal to numel(x)
% 
% INPUTS (optional):
%   s: vector of measure of variability (e.g. standard deviation or
%   confidence interval) for which uncertainty is shown as patch
%   ax: handles of axis for plotting (default = gca)
%   color: color of line and patch (defailt = 'k')
%
% OUTPUTS:
%   line_h: handle to line object
%   patch_h: handle to patch object
%
% Stephen Town - 2016

% Default arguments
if nargin < 5, color = 'k'; end
if nargin < 4, ax = gca; end
if nargin < 3, s = '-'; end

% Hold axes for plotting
axes(ax)
hold on

% Vary approach depending on method
if ~ischar(s)
    
    % Classic approach in which measure of variability is an input
    x2 = [x(:); flipud(x(:))];
    y2 = [y(:)-s(:); flipud(y(:)+s(:))];
    
    patch_h = patch(x2, y2, y2.*0,...   % plot first for layering effect
        'EdgeColor','none',...
        'FaceAlpha',0.5,...
        'FaceColor',color);
    line_h = plot(x,y,'color', color);
    
else
    % User friendly approach where y is provided as a matrix where each row
    % is a sample and mean & standard error across samples are computed
    
    % Calculate mean and standard error
    ym = mean(y,1); % Average across rows
    yn = size(y,1);
    ys = std(y,[],1) ./ sqrt(yn-1);
    
    % Rethrow classic plotting
    [line_h, patch_h] = plotSE_patch(x, ym, ys, ax, color);
end
