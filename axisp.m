function axisp(varargin)
%top_percentage_of_y_for_x_lim limits current x axis based on top y_percentage of y
%   axisp(y_percentage)
%   This functions uses xlim to limit x axis based on the first and the
%   last data points in y that are in the top y_percentage.
%   axisp(y_percentage, 'axis_mode')
%   This function can be used to peform an axis function e.g. axis(tight)
%   before limiting x based on y percentage.
if (nargin < 1) || (nargin > 2)
    keyboard
    error('Wrong number of input arguments')
elseif nargin == 2;
    axis(varargin{2})
end
y_percentage = varargin{1};
h = findobj(gca,'type','line');
y = get(h,'YData');
current_x_lim = get(gca, 'xlim');
y_limit =  max(y).*(1-y_percentage/100) + min(y).*y_percentage/100; % ==  max(y)-y_percentage/100.*(max(y)-min(y))
beyond_y_limit_indecies = find(y > y_limit);
if length(beyond_y_limit_indecies)>1
    future_x_lim = current_x_lim(1)+[beyond_y_limit_indecies(1), beyond_y_limit_indecies(end)]/length(y)*diff(current_x_lim);
    set(gca,'xlim',future_x_lim)
end
end