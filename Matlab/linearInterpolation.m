% Linear interpolation between two values

function [x,y] = linearInterpolation(x_data,y_data,N)

% Initialisation

x(1) = x_data(1);
y(1) = y_data(1);

% Calculation

for i = 1: ( length(x_data) - 1 )
    
    dx = (x_data(i + 1) - x_data(i)) / N;
    dy = (y_data(i + 1) - y_data(i)) / N;
    
    for j = 2 : (N + 1)
        
        x((i - 1) * N + j) = x((i - 1) * N + (j - 1)) + dx;
        y((i - 1) * N + j) = y((i - 1) * N + (j - 1)) + dy;
        
    end
    
end