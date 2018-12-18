function [linehandle, e_handle] = plotAvgWithErrorBar(x_axis, avg, err, freq, color, linewidth)
% plotAvgWithErrorBar
X_indcs =  round((0 : length(x_axis)/freq : length(x_axis)));
% remove the edges
X_indcs = X_indcs( 2 : (end - 1));
X_bar = x_axis(X_indcs);

Y_bar = avg(X_indcs);
Error_bar = err(X_indcs);

% default values
if nargin < 5; linewidth = 1.5; end
if nargin < 4; color = 'b'; end % default color is blue

%% Plot 
linehandle = plot(x_axis, avg, color, 'LineWidth', linewidth); hold on
e_handle = errorbar(X_bar, Y_bar, Error_bar);

e_handle.Marker = 'x';
e_handle.MarkerSize = 5;
e_handle.Color = color;
e_handle.CapSize = 5;
e_handle.LineStyle = 'None';

end
