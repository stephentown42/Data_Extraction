function draw_timeline_heatmap(dt)
%
% Produces a heatmap in the style of the github commit graph in which every
% day is a cell. Here, every year is shown on a separate axis on a separate
% row
%
% Inputs: 
%   - dt: array of datetimes for events
%
% Stephen Town - 7th April 2020

% Get number of years 
years = unique( year( dt));
n_yr = numel( years);

% Create figure and set colormap
figure
cmax = 0;   % For equating colormaps across years
sp = nan( n_yr, 1);

if exist('cmocean','file')
    cmap = cmocean('tempo');
else
    cmap = flipud(colormap('gray'));
end

% For each year
for i = 1 : n_yr 

    sp(i) = subplot(n_yr, 1, i);
    
    draw_year( dt( year(dt) == years(i))); 
    
    cmax = max([cmax get(gca,'clim')]);    
    colormap( sp(i), cmap)
end

% Equate colorscales
set(sp,'clim',[0 cmax])

% Add colorbar
% pass - I'll do this later




function draw_year( dt)
%
% dt: array of datetimes for the specific year

im = zeros(7, 53); 
day_idx = weekday( dt); % Day of the week (1-7)
week_idx = week( dt); % Week of the year (1-53)... 53 = ceil(365/7), don't argue with me!

% For each event (slow but allows for multiple events on one day)
for i = 1 : numel(dt)    
    im(day_idx(i), week_idx(i)) = im(day_idx(i), week_idx(i)) + 1;
end

imagesc(im)

% Make pretty axes and title
set(gca,'xcolor',[0 0 0] + 0.8,...
        'xtick', 0.5:52.5,...
        'xgrid','on',...
        'xticklabel','',...
        'ycolor',[0 0 0] + 0.8,...
        'ygrid','on',...
        'ytick',0.5:6.5,...
        'yticklabel','',...
        'TickLength',[0 0],...
        'gridColor', 'w',...
        'GridAlpha', 1)

ty(1) = text(-1/2, 2, 'Mon');    % Plot as y ticks as text objects
ty(2) = text(-1/2, 4, 'Wed');    % to allow axis color 
ty(3) = text(-1/2, 6, 'Fri');    % to be made invisible
set(ty, 'FontSize',8, 'HorizontalAlignment', 'right')

x_dt = datetime(2020,01,14) + calmonths(0:11);  % Add month x ticks
tx = nan(12,1);
for i = 1 : 12    
    tx(i) = text( week( x_dt(i)), -1/2, datestr(x_dt(i), 'mmm'));
end
set(tx, 'FontSize',8, 'HorizontalAlignment', 'center')

text(-4, -2, sprintf('%d', year(dt(1))),...    % Title
    'FontWeight','bold','FontSize',11)

pbaspect([52 7 1])
